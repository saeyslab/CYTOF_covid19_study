source("R_scripts/load_sample_information.R")
source("R_scripts/load_marker_information.R")
source("R_scripts/plot_files_scatters.R")
source("R_scripts/plot_2D_scatters.R")
source("R_scripts/make_aggregate.R")
source("R_scripts/quantiles.R")
source("R_scripts/dendrogram.R")
library(FlowSOM)
library(ggrepel)

dir.create("Results")
dir.create("RDS")

# Load sample information ------------------------------------------------------

recompute <- FALSE

samples_meta_file <- "Meta_information/Extracted_Data_24042020.xlsx"
selection_file <- "Meta_information/Template Contagious Data_CyTOF.xlsx"
PeacoQC_fcs_dir <- "PeacoQC_Data_05_05/PeacoQC_results/fcs_files"

samples <- load_sample_information(meta_file = samples_meta_file,
                                   selection_file = selection_file,
                                   fcs_dir = PeacoQC_fcs_dir,
                                   recompute = recompute)

samples <- samples[order(as.POSIXct(samples$Date_measured, 
                                    format = "%d-%b-%Y %H:%M:%S")), ]

plotSampleOverview(samples,
                   plot = "Results/SampleOverview.pdf")

# Samples to use in this study -------------------------------------------------

## A - For training ------------------------------------------------------------

#test <- read.FCS("Aggregate/aggregate_selectedForFlowSOM.fcs")


passed_peacoQC <- !is.na(samples$File_peacoQC)
enough_cells <- samples$nCells >= 50000
only_dip_panel <- c("COVID_ICU_013 4_ICU_A",
                    "COVID_ICU_014 4_ICU_A",
                    "COVID_ICU_015 4_ICU_A",
                    "COVID_ICU_016 4_ICU_A")
pca_outliers <- c("COVID_ICU_003 8_ICU_D",
                  "COVID_W_020 1_W_O",
                  "COVID_W_023 1_W_O",
                  "COVID_HC_002 0_CTR")
controls <- grep("HEALTHY", samples$Sample_ID, value = TRUE)
doubles <-  grep("\\.1", samples$Sample_ID, value = TRUE)
clinical_outliers <- c("COVID_W_005 1_W_O",
                       "COVID_W_013 1_W_O",
                       "COVID_ICU_048 4_ICU_A")

use_to_train <- passed_peacoQC & enough_cells & 
  ! samples$Sample_ID %in% c(only_dip_panel, pca_outliers,
                             controls, doubles)

samples_for_training <- samples[use_to_train, ]

## B - For statistics ----------------------------------------------------------
samples[pca_outliers, "QC_cytof"] <- "0"

samples_of_interest <- samples[
  intersect(
    intersect(
      intersect(c(grep("healthy", samples$Clinical_condition_at_sampling),
                  grep("good", samples$Clinical_condition_at_sampling),
                  grep("bad", samples$Clinical_condition_at_sampling)),
                grep("1", samples$`measured_in_cytof?`)),
      grep("1", samples$QC_cytof)),
    which(! samples$Sample_ID %in% clinical_outliers)),]

# Markers to use in the study --------------------------------------------------

markers_measured <- load_markers_measured(samples,
                                          excel_file = "Meta_information/markers_immune_cytof_panel_v4 (1).xlsx",
                                          sheet_name = "Extended panel",
                                          recompute = recompute)

markers_df <- dplyr::filter(markers_measured,
                            markers_measured$channel_type == "DIP panel" |
                              markers_measured$Marker %in% c("162Dy_CD69",
                                                             "165Ho_CD163",
                                                             "169Tm_NKG2A"))
markers <- unique(markers_df$Marker)
names(markers) <- unique(markers_df$Channel)
markers <- markers[order(as.numeric(gsub("..([0-9]*)Di", "\\1", names(markers))))]
length(markers)

marker_info <- load_marker_excel()
DIP_channels  <- marker_info$Channel[grep("DIP", marker_info$`CyTOF panel`)]
phenotypic_channels  <- marker_info$Channel[grep("phenotypic", marker_info$`marker type`)]
phenotypic_channels <- phenotypic_channels[-c(1:3)]

# FlowSOM Step 1 ---------------------------------------------------------------
dir.create("Aggregate")

aggregate_step1 <- make_aggregate(samples = samples_for_training,
                                    channels = names(markers),
                                    cTotal = 3000000,
                                    agg_file = "Aggregate/aggregate_step1.fcs",
                                    seed = 1,
                                    recompute = recompute)

marker_subset <- sapply(c("CD45$", "CD66b$", "CD3$", "CD4$",
                          "CD8a$", "TCRgd$", "CD56_NCAM$", "CD11c$",
                          "CD19$", "CD14$", "CD20$"),
                        grep, markers, value = TRUE)
names(marker_subset) <- gsub(".*\\.", "", names(marker_subset))

n_meta <- 30

fsom_file <- "RDS/FlowSOM_step1.RDS"
if(!recompute & file.exists(fsom_file)){
  fsom <- readRDS(fsom_file)
  message("Reloaded FlowSOM result. Set recompute = TRUE to recompute.")
} else {
  fsom <- FlowSOM(aggregate_step1,
                  colsToUse = names(marker_subset),
                  scale = FALSE,
                  nClus = n_meta,
                  seed = 1)
  fsom$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", "\\1", fsom$FlowSOM$prettyColnames[1:33])
  saveRDS(fsom, fsom_file)
}

metacluster_mfis <- MetaclusterMFIs(fsom)
metacluster_mfis <- metacluster_mfis[, names(markers)]
colnames(metacluster_mfis) <- markers
rownames(metacluster_mfis) <-  paste0("MC",rownames(metacluster_mfis))

selected_mc <- which(metacluster_mfis[,"172Yb_CD66b"] < 2 | metacluster_mfis[,"89Y_CD45"] > 4 )



dir.create("NonGranulocyte")
dir.create("Granulocyte")
dir.create("Results/splitGranulocytes")
for (sample in c(rownames(samples_for_training), only_dip_panel)) {

  message("Mapping ",sample)
  if (sample %in% only_dip_panel){
    ff <- read.FCS(samples_of_interest[sample, "File_peacoQC"])
  } else {
    ff <- read.FCS(samples_for_training[sample, "File_peacoQC"])
  }
  
  fsom_tmp <- NewData(fsom, ff)
  selection <- GetMetaclusters(fsom_tmp) %in% selected_mc
  
  if (sample %in% only_dip_panel){
    write.FCS(ff[selection,],
              file.path("NonGranulocyte", 
                        basename(samples_of_interest[sample, "File_peacoQC"])))
    
    write.FCS(ff[!selection,],
              file.path("Granulocyte", 
                        basename(samples_of_interest[sample, "File_peacoQC"])))
  } else {
    write.FCS(ff[selection,],
              file.path("NonGranulocyte", 
                        basename(samples_for_training[sample, "File_peacoQC"])))
    
    write.FCS(ff[!selection,],
              file.path("Granulocyte", 
                        basename(samples_for_training[sample, "File_peacoQC"])))
  }
  
  metacluster_mfis_tmp <- MetaclusterMFIs(fsom_tmp)
  colnames(metacluster_mfis_tmp) <- get_markers(ff, colnames(metacluster_mfis_tmp))
  rownames(metacluster_mfis_tmp) <-  paste0("MC",rownames(metacluster_mfis_tmp))
  
  png(paste0("Results/splitGranulocytes/",sample,".png"),
      width = 1500, height = 1000)
  plot(ff@exprs[1:20000, 
                get_channels(ff, c("CD66b", "CD45"))],
       pch = ".",
       col = ifelse(GetMetaclusters(fsom_tmp)[1:20000] %in% selected_mc,
                    "#7bccc4",
                    "black"),
       cex = 10, cex.axis = 2, cex.lab = 2, cex.main = 2,
       xlab = "CD66b", ylab = "CD45",
       xlim = c(0,8), ylim = c(0,8))
  points(metacluster_mfis_tmp[,c("172Yb_CD66b", "89Y_CD45")],
         pch = 19,
         col = ifelse(seq_len(n_meta) %in% selected_mc, "#2b8cbe", "#fb6a4a"),
         cex =  3)
  abline(v = 2, col = "grey", lwd = 3)
  abline(h = 4, col = "grey", lwd = 3)
  dev.off()
}

samples_for_training$File_peacoQC <- file.path("NonGranulocyte",
                                               basename(samples_for_training$File_peacoQC))

samples_of_interest$File_peacoQC <- file.path("NonGranulocyte",
                                              basename(samples_of_interest$File_peacoQC))

# FlowSOM Step 2 ---------------------------------------------------------------

aggregate_nonGranulocyte <- make_aggregate(samples = samples_for_training,
                                           channels = names(markers),
                                           cTotal = 3000000,
                                           agg_file = "Aggregate3/aggregate_step2.fcs",
                                           seed = 1,
                                           recompute = recompute)

fsom_step2_file <- "RDS/FlowSOM_step2_repeat2.RDS"
if(!recompute & file.exists(fsom_step2_file)){
  fsom_step2 <- readRDS(fsom_step2_file)
  message("Reloaded FlowSOM result. Set recompute = TRUE to recompute.")
} else {
  fsom_step2 <- FlowSOM(aggregate_nonGranulocyte,
                        colsToUse = names(markers),
                        scale = FALSE,
                        nClus = n_meta,
                        seed = 1)
  fsom_step2$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", 
                                                  "\\1", 
                                                  fsom_step2$FlowSOM$prettyColnames[1:33])
  saveRDS(fsom_step2, fsom_step2_file)
}


selection_cl30 <- GetClusters(fsom_step2) == 30
ggplot(data.frame(aggregate_nonGranulocyte@exprs[selection_cl30,])) +
  geom_point(aes(x=Pr141Di, 
                 y=Eu153Di, 
                 color=samples_for_training$Sample_ID[File]),
             size = 0.5) +
  theme_minimal() +
  xlim(0,8) + ylim(0,8) +
  theme(legend.position = "none")

sort(counts[,30],decreasing = TRUE)

fsom_step2_extended_file <- "RDS/FlowSOM_step2_repeat2_extended.RDS"
if(!recompute & file.exists(fsom_step2_extended_file)){
  fsom_step2 <- readRDS(fsom_step2_extended_file)
  message("Reloaded previous result, put recompute to TRUE to recompute.")
} else {
  manually_adapted_clusters <- list(list(original = 30,
                                         "Eu153Di" = c(2, 5)),
                                    list(original = 14,
                                         "Eu153Di" = c(2, 5),
                                         "Pr141Di" = c(5, 7)),
                                    list(original = 15,
                                         "Gd156Di" = c(4, 4, 0),
                                         "Sm152Di" = c(3, 1, 3)),
                                    list(original = 25,
                                         "Gd156Di" = c(4, 4),
                                         "Sm152Di" = c(3, 1)),
                                    list(original = 26,
                                         "Gd156Di" = c(4, 4),
                                         "Sm152Di" = c(3, 1)),
                                    list(original = 7,
                                         "Eu151Di" = c(2.5, 0)))
  for(cl_id in seq_along(manually_adapted_clusters)){
    cl_info <- manually_adapted_clusters[[cl_id]]
    original_cluster <- cl_info[["original"]]
    channels <- names(cl_info[-1])
    
    for(channel in channels){
      fsom_step2$FlowSOM$map$codes[original_cluster, channel] <- cl_info[[channel]][1]
    }
    
    for(i in seq_len(length(cl_info[[channel]])-1)){
      copy <- fsom_step2$FlowSOM$map$codes[original_cluster,]
      for(channel in channels){
        copy[channel] <- cl_info[[channel]][i+1]
      }
      fsom_step2$FlowSOM$map$codes <- rbind(fsom_step2$FlowSOM$map$codes,
                                            copy)
      fsom_step2$FlowSOM$map$nNodes <- fsom_step2$FlowSOM$map$nNodes + 1
      
    }
    
  }
  
  rownames(fsom_step2$FlowSOM$map$codes) <- NULL
  
  fsom_step2$FlowSOM$map$grid <- rbind(fsom_step2$FlowSOM$map$grid,
                                       cbind("Var1" = 11, "Var2" = 1:7))
  
  fsom_step2$FlowSOM <- FlowSOM::BuildMST(fsom_step2$FlowSOM)
  fsom_step2 <- NewData(fsom_step2, aggregate_nonGranulocyte)
  fsom_step2$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", 
                                                  "\\1", 
                                                  fsom_step2$FlowSOM$prettyColnames[1:33])
  
  saveRDS(fsom_step2, fsom_step2_extended_file)
}

fsom_step2_sub <- fsom_step2
fsom_step2_sub$FlowSOM$map$codes <- fsom_step2_sub$FlowSOM$map$codes[, DIP_channels]


# FlowSOM Step 2 Granulocyte ---------------------------------------------------


samples_for_training$File_peacoQC <- file.path("Granulocyte",
                                               basename(samples_for_training$File_peacoQC))

samples_of_interest$File_peacoQC <- file.path("Granulocyte",
                                              basename(samples_of_interest$File_peacoQC))

aggregate_granulocyte <- make_aggregate(samples = samples_for_training,
                                        channels = names(markers),
                                        cTotal = 3000000,
                                        agg_file = "Aggregate/aggregate_step2_granulocyte.fcs",
                                        seed = 1,
                                        recompute = recompute)

fsom_step2_file_gr <- "RDS/FlowSOM_step2_gr.RDS"
if(!recompute & file.exists(fsom_step2_file_gr)){
  fsom_step2_gr <- readRDS(fsom_step2_file_gr)
  message("Reloaded FlowSOM result. Set recompute = TRUE to recompute.")
} else {
  fsom_step2_gr <- FlowSOM(aggregate_granulocyte,
                           colsToUse = names(markers),
                           scale = FALSE,
                           nClus = n_meta,
                           seed = 1)
  fsom_step2_gr$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", 
                                                  "\\1", 
                                                  fsom_step2_gr$FlowSOM$prettyColnames[1:33])
  saveRDS(fsom_step2_gr, fsom_step2_file_gr)
}

PlotFlowSOM(fsom_step2_gr)

fsom_step2_sub_gr <- fsom_step2_gr
fsom_step2_sub_gr$FlowSOM$map$codes <- fsom_step2_sub_gr$FlowSOM$map$codes[, DIP_channels]


# Mapping samples of interest --------------------------------------------------

counts_file <- "RDS/counts.RDS"
mfis_file <- "RDS/mfis.RDS"

if(file.exists(counts_file)){
  counts <- readRDS(counts_file)
  mfis <- readRDS(mfis_file)
} else {
  counts <- matrix(0,
                   nrow = length(samples_of_interest$Sample_ID),
                   ncol = fsom_step2$FlowSOM$map$nNodes + 
                     fsom_step2_gr$FlowSOM$map$nNodes,
                   dimnames = list(samples_of_interest$Sample_ID,
                                   c(paste0("Cl", seq_len(fsom_step2$FlowSOM$map$nNodes)),
                                     paste0("grCl", seq_len(fsom_step2_gr$FlowSOM$map$nNodes)))))
  mfi_names <- expand.grid(c(paste0("Cl", seq_len(fsom_step2$FlowSOM$map$nNodes)),
                             paste0("grCl", seq_len(fsom_step2_gr$FlowSOM$map$nNodes))),
                           markers)
  mfi_names <- paste0(mfi_names[,1], " ", mfi_names[,2])
  mfis <- matrix(NA,
                 nrow = length(samples_of_interest$Sample_ID),
                 ncol = length(mfi_names),
                 dimnames = list(samples_of_interest$Sample_ID,
                                 mfi_names))
  plots <- list()
  for (sample in rownames(samples_of_interest)) {
    message("Mapping ",sample)
    
    # Non Granulocyte --- 
    
    ff <- read.FCS(file.path("NonGranulocyte",
                             basename(samples_of_interest[sample, "File_peacoQC"])))
    if (sample %in% only_dip_panel){
      fsom_tmp <- NewData(fsom_step2_sub, ff)
    } else {
      fsom_tmp <- NewData(fsom_step2, ff)
    }
    fsom_tmp$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", "\\1", fsom_tmp$FlowSOM$prettyColnames[1:33])
    
    if(!sample %in% only_dip_panel){
      plots[[sample]] <- PlotFlowSOM(fsom_tmp$FlowSOM,
                                    title = sample)
    } else {
      plots[[sample]] <- PlotFlowSOM(fsom_tmp,
                                     markers = DIP_channels,
                                     title = sample)
    }
    
    counts_tmp <- table(GetClusters(fsom_tmp))
    counts[sample, paste0("Cl", names(counts_tmp))] <- counts_tmp
    
    if (sample %in% only_dip_panel){
      mfis_tmp <- data.frame(Cluster = paste0("Cl", seq_len(fsom_step2$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_tmp)[,DIP_channels])
      colnames(mfis_tmp) <- c("Cluster", markers[DIP_channels])
    } else {
      mfis_tmp <- data.frame(Cluster = paste0("Cl", seq_len(fsom_step2$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_tmp)[,names(markers)])
      colnames(mfis_tmp) <- c("Cluster", markers)
    }
    mfis_tmp <- gather(mfis_tmp, "Marker", "MFI", -Cluster)
    mfis[sample, paste0(mfis_tmp$Cluster," ",mfis_tmp$Marker)] <- mfis_tmp$MFI
    
    # Granulocyte --- 
    
    ff <- read.FCS(file.path("Granulocyte",
                             basename(samples_of_interest[sample, "File_peacoQC"])))
    if (sample %in% only_dip_panel){
      fsom_tmp <- NewData(fsom_step2_sub_gr, ff)
    } else {
      fsom_tmp <- NewData(fsom_step2_gr, ff)
    }
    fsom_tmp$FlowSOM$prettyColnames[1:33] <- gsub(".*_(.*) <.*", "\\1", fsom_tmp$FlowSOM$prettyColnames[1:33])
    
    if(!sample %in% only_dip_panel){
      plots[[paste0(sample, "_gr")]] <- PlotFlowSOM(fsom_tmp,
                                                    title = sample)
    } else {
      plots[[paste0(sample, "_gr")]] <- PlotFlowSOM(fsom_tmp,
                                                    markers = DIP_channels,
                                                    title = sample)
    }
    
    counts_tmp <- table(GetClusters(fsom_tmp))
    counts[sample, paste0("grCl", names(counts_tmp))] <- counts_tmp
    
    if (sample %in% only_dip_panel){
      mfis_tmp <- data.frame(Cluster = paste0("grCl", seq_len(fsom_step2_gr$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_tmp)[,DIP_channels])
      colnames(mfis_tmp) <- c("Cluster", markers[DIP_channels])
    } else {
      mfis_tmp <- data.frame(Cluster = paste0("grCl", seq_len(fsom_step2_gr$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_tmp)[,names(markers)])
      colnames(mfis_tmp) <- c("Cluster", markers)
    }
    mfis_tmp <- gather(mfis_tmp, "Marker", "MFI", -Cluster)
    mfis[sample, paste0(mfis_tmp$Cluster," ",mfis_tmp$Marker)] <- mfis_tmp$MFI
  }
  
  pdf("Results/individual_samples.pdf",
      width = 11.7,
      height = 8.3)
  invisible(lapply(plots, plot))
  dev.off()
  saveRDS(counts, counts_file)
  saveRDS(mfis, mfis_file)
}


# Derived values ---------------------------------------------------------------

pctgs <- counts
pctgs[,grep("^Cl", colnames(pctgs))] <- t(apply(counts[,grep("^Cl", colnames(pctgs))], 1, function(x) x / sum(x) * 100))
pctgs[,grep("^grCl", colnames(pctgs))] <- t(apply(counts[,grep("^grCl", colnames(pctgs))], 1, function(x) x / sum(x) * 100))

hierarchical_clustering <- hclust(dist(fsom_step2$FlowSOM$map$codes))
hierarchical_clustering$labels <- paste0("Cl", 1:107)
merge_groups <- list()
for (merge_id in 1:nrow(hierarchical_clustering$merge)){
  merge_groups[[merge_id]] <- list()
  groups <- hierarchical_clustering$merge[merge_id,]
  for(side in 1:2){
    if(groups[side] < 0){
      merge_groups[[merge_id]][[side]] <- -groups[side]
    } else {
      merge_groups[[merge_id]][[side]] <- c(merge_groups[[groups[side]]][[1]], 
                                            merge_groups[[groups[side]]][[2]])
    }
  }
}

ratios <- matrix(NA,
                 nrow = nrow(pctgs),
                 ncol = length(merge_groups),
                 dimnames = list(rownames(pctgs), 
                                 seq_len(length(merge_groups))))
parent_clusters <- matrix(NA,
                          nrow = nrow(pctgs),
                          ncol = length(merge_groups),
                          dimnames = list(rownames(pctgs), 
                                          seq_len(length(merge_groups))))
for(merge_id in seq_along(merge_groups)){
  colnames(ratios)[merge_id] <- 
    paste("(", 
          paste(hierarchical_clustering$labels[merge_groups[[merge_id]][[1]]], 
                collapse = " + "),
          ") / (",
          paste(hierarchical_clustering$labels[merge_groups[[merge_id]][[2]]],
                collapse = " + "), 
          ")")
  colnames(parent_clusters)[merge_id] <- 
    paste(c(hierarchical_clustering$labels[merge_groups[[merge_id]][[1]]],
            hierarchical_clustering$labels[merge_groups[[merge_id]][[2]]]),
          collapse = " + ")
  side1 <- rowSums(pctgs[,merge_groups[[merge_id]][[1]], drop = FALSE]) + 1e-100
  side2 <- rowSums(pctgs[,merge_groups[[merge_id]][[2]], drop = FALSE]) + 1e-100

  parent_clusters[,merge_id] <- side1+side2
  ratios[,merge_id] <- side1 / side2
}


hierarchical_clustering_gr <- hclust(dist(fsom_step2_gr$FlowSOM$map$codes))
hierarchical_clustering_gr$labels <- paste0("grCl", 1:100)
merge_groups_gr <- list()
for (merge_id in 1:nrow(hierarchical_clustering_gr$merge)){
  merge_groups_gr[[merge_id]] <- list()
  groups <- hierarchical_clustering_gr$merge[merge_id,]
  for(side in 1:2){
    if(groups[side] < 0){
      merge_groups_gr[[merge_id]][[side]] <- -groups[side]
    } else {
      merge_groups_gr[[merge_id]][[side]] <- c(merge_groups_gr[[groups[side]]][[1]], 
                                               merge_groups_gr[[groups[side]]][[2]])
    }
  }
}

ratios_gr <- matrix(NA,
                 nrow = nrow(pctgs),
                 ncol = length(merge_groups_gr),
                 dimnames = list(rownames(pctgs), 
                                 seq_len(length(merge_groups_gr))))
parent_clusters_gr <- matrix(NA,
                          nrow = nrow(pctgs),
                          ncol = length(merge_groups_gr),
                          dimnames = list(rownames(pctgs), 
                                          seq_len(length(merge_groups_gr))))
for(merge_id in seq_along(merge_groups_gr)){
  colnames(ratios_gr)[merge_id] <- 
    paste("(", 
          paste(hierarchical_clustering_gr$labels[merge_groups_gr[[merge_id]][[1]]], 
                collapse = " + "),
          ") / (",
          paste(hierarchical_clustering_gr$labels[merge_groups_gr[[merge_id]][[2]]],
                collapse = " + "), 
          ")")
  colnames(parent_clusters_gr)[merge_id] <- 
    paste(c(hierarchical_clustering_gr$labels[merge_groups_gr[[merge_id]][[1]]],
            hierarchical_clustering_gr$labels[merge_groups_gr[[merge_id]][[2]]]),
          collapse = " + ")
  side1 <- rowSums(pctgs[,merge_groups_gr[[merge_id]][[1]], drop = FALSE]) + 1e-100
  side2 <- rowSums(pctgs[,merge_groups_gr[[merge_id]][[2]], drop = FALSE]) + 1e-100
  
  parent_clusters_gr[,merge_id] <- side1+side2
  ratios_gr[,merge_id] <- side1 / side2
}

# Save to excel ----------------------------------------------------------------

xlsx::write.xlsx(counts,
                 file = "Results/FlowSOM_results.xlsx",
                 sheetName = "Counts")

xlsx::write.xlsx(pctgs,
                 file = "Results/FlowSOM_results.xlsx",
                 sheetName = "Percentages",
                 append = TRUE)

xlsx::write.xlsx(data.frame(ratios, 
                            ratios_gr,
                            check.names = FALSE),
                 file = "Results/FlowSOM_results.xlsx",
                 sheetName = "Ratios",
                 append = TRUE)

xlsx::write.xlsx(data.frame(parent_clusters, 
                            parent_clusters_gr,
                            check.names = FALSE),
                 file = "Results/FlowSOM_results.xlsx",
                 sheetName = "Parent clusters",
                 append = TRUE)

xlsx::write.xlsx(data.frame(Total = rowSums(counts), 
                            NonGranulocytes = rowSums(counts[,grep("^Cl", colnames(counts))]),
                            Granulocytes = rowSums(counts[,grep("^grCl", colnames(counts))]),
                            check.names = FALSE),
                 file = "Results/FlowSOM_results.xlsx",
                 sheetName = "Total counts",
                 append = TRUE)

xlsx::write.xlsx(data.frame(t(mfis),
                            check.names = FALSE),
                 file = "Results/FlowSOM_mfis.xlsx",
                 sheetName = "MFIs")

# Cluster colors ---------------------------------------------------------------

cluster_annotation <- xlsx::read.xlsx("Meta_information/FlowSOM_cluster_annotations_v2.xlsx",
                                      sheetName = "Sheet2",
                                      check.names = FALSE)
rownames(cluster_annotation) <- cluster_annotation$cluster
cluster_label <- cluster_annotation[c(paste0("Cl", 1:107), paste0("grCl", 1:100)), "higher level"]

cluster_label_final <- cluster_annotation[c(paste0("Cl", 1:107), paste0("grCl", 1:100)), "Final_annotation"]

manual_colors <- c("mDC" = '#fff79e',
                   "monocytes" = '#88419d',
                   "CD4 T cells" = '#dd3497',
                   "T cells DN" = '#31a354',
                   "CD8 T cells" = '#66c2a4',
                   "neutrophil" = "yellow",
                   "CD8 NKT" = "#ee204d",
                   "NK cells" = '#000080',
                   "undefined" = 'grey',
                   "gdTcells" = '#ccff00',
                   "B cells" = '#a6bddb',
                   "Basophils" = '#993404',
                   "pDC" = '#d492e8',
                   "eosinophils" = '#ffcc99')

cluster_color <- manual_colors[cluster_label]

hierarchical_clustering$labels <- cluster_annotation[paste0("Cl",1:107), "Final_annotation"]
hierarchical_clustering_gr$labels <- cluster_annotation[paste0("grCl",1:100), "Final_annotation"]

# Figures ----------------------------------------------------------------------

p <- PlotVariable(fsom_step2,
                  equalNodeSize = TRUE,
                  maxNodeSize = 0.25,
                  factor(cluster_label[1:107], levels = names(manual_colors)),
                  colorPalette = manual_colors)
layout <- data.frame(fsom_step2$FlowSOM$MST$l)
colnames(layout) <- c("V1", "V2")
p <-  AddLabels(p, 
                layout = layout, 
                metadata = 1:107)
p
ggsave("Results/tree_overview.pdf", p,
       width = 10, height = 10)

p <- PlotVariable(fsom_step2_gr,
                  equalNodeSize = TRUE,
                  maxNodeSize = 0.25,
                  factor(cluster_label[108:207], levels = names(manual_colors)),
                  colorPalette = manual_colors)
layout <- data.frame(fsom_step2_gr$FlowSOM$MST$l)
colnames(layout) <- c("V1", "V2")
p <-  AddLabels(p, 
                layout = layout, 
                metadata = 1:100)
p

ggsave("Results/tree_overview_gr.pdf", p,
       width = 10, height = 10)

cluster_mfis <- GetMFIs(fsom_step2, colsUsed = TRUE, prettyColnames = TRUE)
cluster_mfis_gr <- GetMFIs(fsom_step2_gr, colsUsed = TRUE, prettyColnames = TRUE)

column_order <- c('CD45',  
                  'CD66b',
                  'CD11c',  
                  'HLA-DR',
                  'CD19', 'CD20', 'IgD', 
                  'CD14', 'CD16', 
                  'CD3',  'TCRgd', 'CD4', 'CD8a', 'CCR4','CD45RO', 'CD45RA',
                  'NCAM',  'CD161',
                  'CD27',  'CD28', 'CD38',  'CD57', 'CD69',  'CD163', 'CD294',
                  'CCR6', 'CCR7', 'CXCR3', 'CXCR5',  'IL-2Ra', 'IL-3R', 'IL-7Ra', 'NKG2A')

heatmap_cl_mfi <- pheatmap::pheatmap(cluster_mfis[hierarchical_clustering$order,
                                                  column_order],
                                     cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     labels_row = hierarchical_clustering$labels[hierarchical_clustering$order],
                                     main = "Cluster MFIs",
                                     annotation_row = data.frame(annotation = factor(cluster_label, 
                                                                                     levels = names(manual_colors))),
                                     annotation_colors = list(annotation = manual_colors),
                                     annotation_legend = FALSE)

ggsave(filename = "Results/cluster_heatmap.pdf",
       plot = heatmap_cl_mfi,
       width = 15, height = 15)

heatmap_cl_mfi_gr <- pheatmap::pheatmap(cluster_mfis_gr[hierarchical_clustering_gr$order,
                                                        column_order],
                                     cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     labels_row = hierarchical_clustering_gr$labels[hierarchical_clustering_gr$order],
                                     main = "Granulocyte Cluster MFIs",
                                     annotation_row = data.frame(annotation = factor(cluster_label[108:207], 
                                                                                     levels = names(manual_colors))),
                                     annotation_colors = list(annotation = manual_colors),
                                     annotation_legend = FALSE)

ggsave(filename = "Results/cluster_heatmap_gr.pdf",
       plot = heatmap_cl_mfi_gr,
       width = 15, height = 15)

plots <- list()
for(channel in sapply(column_order, function(m) names(markers)[grep(paste0(m,"$"), markers)])){
  plots[[channel]] <- PlotFlowSOM(fsom_step2,
                                  equalNodeSize = TRUE,
                                  maxNodeSize = 0.25,
                                  type = "marker",
                                  metadata = channel,
                                  plot = FALSE) + 
    ggtitle(markers[[channel]]) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8))
}

ggsave("Results/FlowSOM_per_marker.pdf",
       plot = do.call(ggarrange, c(plots, list(common.legend = TRUE))),
       width = 11.7, height = 8.3)

plots <- list()
for(channel in sapply(column_order, function(m) names(markers)[grep(paste0(m,"$"), markers)])){
  plots[[channel]] <- PlotFlowSOM(fsom_step2_gr,
                                  equalNodeSize = TRUE,
                                  maxNodeSize = 0.25,
                                  type = "marker",
                                  metadata = channel,
                                  plot = FALSE) + 
    ggtitle(markers[[channel]]) +
    theme(legend.position = "none",
          plot.title = element_text(size = 8))
}

ggsave("Results/FlowSOM_per_marker_gr.pdf",
       plot = do.call(ggarrange, c(plots, list(common.legend = TRUE))),
       width = 11.7, height = 8.3)

# Volcano plot -----------------------------------------------------------------

all_features <- cbind(pctgs, ratios, ratios_gr, parent_clusters, parent_clusters_gr, mfis)

volcano_plots <- list()

for(groups in list(c("good", "bad"),
                   c("healthy", "good"),
                   c("healthy", "bad"))){
  group1 <- groups[1]
  group2 <- groups[2]
  group_name <- paste0(group1, " vs ", group2)

  p_values <- apply(all_features, 2, function(feature){
    values_group1 <- feature[samples_of_interest$Clinical_condition_at_sampling == group1]
    values_group2 <- feature[samples_of_interest$Clinical_condition_at_sampling == group2]
    if (! (all(is.na(values_group1)) | all(is.na(values_group2)))) {
      test_res <- wilcox.test(values_group1, values_group2, exact = FALSE)
      test_res$p.value
    } else {
      1
    }
    
  })
  
  fold_changes <- apply(all_features, 2, function(feature){
    medians <- c(median(feature[samples_of_interest$Clinical_condition_at_sampling == group1], na.rm = TRUE) + 1e-100,
                 median(feature[samples_of_interest$Clinical_condition_at_sampling == group2], na.rm = TRUE) + 1e-100)
    fold_change <- max(medians) / min(medians) * (-1)^which.max(medians)
      
  })
  
  to_plot <- data.frame("feature" = colnames(all_features),
                        "pvalue" = p_values,
                        "-log10p" = -log10(p_values),
                        "foldchange" = fold_changes,
                        "log10foldchange" = sign(fold_changes)*log10(abs(fold_changes)), #Double check
                        "color" = c(cluster_color,
                                    rep("black", ncol(ratios)),
                                    rep("black", ncol(ratios_gr)),
                                    rep("grey", ncol(parent_clusters)),
                                    rep("grey", ncol(parent_clusters_gr)),
                                    rep("gold", ncol(mfis))),
                        "group" = group_name,
                        "subset" = c(rep("cluster", ncol(pctgs)),
                                     rep("ratio", ncol(ratios)),
                                     rep("ratio", ncol(ratios_gr)),
                                     rep("sum", ncol(parent_clusters)),
                                     rep("sum", ncol(parent_clusters_gr)),
                                     rep("mfi", ncol(mfis))),
                        check.names = FALSE)

  saveRDS(to_plot,
          file = paste0("RDS/p_values_",group_name,".RDS"))
  xlsx::write.xlsx(to_plot[head(order(to_plot$`-log10p`, decreasing=TRUE), n = 100), ],
                   file = paste0("Results/p_values_",group_name,".xlsx"))
  
  volcano_plots[[group_name]] <- ggplot(to_plot[rev(seq_len(nrow(to_plot))),], 
                                        aes(x = log10foldchange, y = `-log10p`)) +
    geom_point(aes(col = color)) +
    geom_text_repel(aes(label = ifelse(`-log10p` > 3 | log10foldchange > 3.5 | log10foldchange < -3.5, feature, "")), 
                    size = 2, min.segment.length = 0) +
    scale_color_identity() +
    xlim(-7, 7) +
    ylim(0, 5) +
    ggtitle(group_name) +
    theme_minimal() +
    facet_grid(group ~ subset)
  
  volcano_plots[[paste0(group_name, "2")]] <- ggplot(to_plot[rev(seq_len(nrow(to_plot))),], 
                                        aes(x = log10foldchange, y = `-log10p`)) +
    geom_point(aes(col = color)) +
    scale_color_identity() +
    xlim(-7, 7) +
    ylim(0, 5) +
    ggtitle(group_name) +
    theme_minimal() +
    facet_grid(group ~ subset)
  
  pdf(paste0("Results/dendrograms ",group_name,".pdf"),
      height = 10,
      width = 10)
  par(mar = c(2,4,2,10))
  
  dendrogram_values <- list()
  dendrogram_values[[paste0(group_name, " ratio", " -log10(p)")]] <- -log10(p_values[to_plot$subset == "ratio"][1:106])
  dendrogram_values[[paste0(group_name, " sum", " -log10(p)")]] <- -log10(p_values[to_plot$subset == "sum"][1:106])
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] <- sign(fold_changes[to_plot$subset == "ratio"][1:106]) *
                                                                           log10(abs(fold_changes[to_plot$subset == "ratio"][1:106]))
  dendro_outlier <- dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] > 3
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]][dendro_outlier] <- 3
  dendro_outlier <- dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] < -3
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]][dendro_outlier] <- -3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] <- sign(fold_changes[to_plot$subset == "sum"][1:106]) *
                                                                         log10(abs(fold_changes[to_plot$subset == "sum"][1:106]))
  dendro_outlier <- dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] > 3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]][dendro_outlier] <- 3
  dendro_outlier <- dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] < -3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]][dendro_outlier] <- -3
  
  dendrogram_values[[]]
  for(dendrogram_group in names(dendrogram_values)){
    plot_cluster_dendrogram(hierarchical_clustering,
                            values = dendrogram_values[[dendrogram_group]],
                            values_title = "",
                            leaf_col = cluster_color[1:107],
                            main = dendrogram_group)
  }
  
  dendrogram_values <- list()
  dendrogram_values[[paste0(group_name, " ratio", " -log10(p)")]] <- -log10(p_values[to_plot$subset == "ratio"][107:205])
  dendrogram_values[[paste0(group_name, " sum", " -log10(p)")]] <- -log10(p_values[to_plot$subset == "sum"][107:205])
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] <- sign(fold_changes[to_plot$subset == "ratio"][107:205]) *
    log10(abs(fold_changes[to_plot$subset == "ratio"][107:205]))
  dendro_outlier <- dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] > 3
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]][dendro_outlier] <- 3
  dendro_outlier <- dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]] < -3
  dendrogram_values[[paste0(group_name, " ratio", " log10foldchange")]][dendro_outlier] <- -3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] <- sign(fold_changes[to_plot$subset == "sum"][107:205]) *
    log10(abs(fold_changes[to_plot$subset == "sum"][107:205]))
  dendro_outlier <- dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] > 3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]][dendro_outlier] <- 3
  dendro_outlier <- dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]] < -3
  dendrogram_values[[paste0(group_name, " sum", " log10foldchange")]][dendro_outlier] <- -3
  
  for(dendrogram_group in names(dendrogram_values)){
    plot_cluster_dendrogram(hierarchical_clustering_gr,
                            values = dendrogram_values[[dendrogram_group]],
                            values_title = "",
                            leaf_col = cluster_color[108:207],
                            main = dendrogram_group)
  }
  dev.off()

  group1_ids <- rownames(samples_of_interest)[samples_of_interest$Clinical_condition_at_sampling == group1]
  group2_ids <- rownames(samples_of_interest)[samples_of_interest$Clinical_condition_at_sampling == group2]
  
  features_top100 <- names(head(sort(p_values), n = 100))
  boxplots <- list()
  for(feature in features_top100){
    boxplots[[feature]] <- ggplot(data.frame(Patient = gsub(" .*", "", c(group1_ids, group2_ids)),
                                             Group = factor(c(rep(group1, length(group1_ids)),
                                                              rep(group2, length(group2_ids))),
                                                            levels = c(group1, group2)),
                                             Value = all_features[c(group1_ids, group2_ids), feature])) +
      geom_boxplot(aes(x = Group, y = Value)) +
      geom_text(aes(x = Group, y = Value, label = Patient), size = 1.5,
                position = ggbeeswarm::position_quasirandom()) +
      ggtitle(paste0(feature, ": p ", format(p_values[feature], digits = 2),
                     " fc ", format(fold_changes[feature], digits = 2))) +
      theme_minimal() +
      theme(title = element_text(size = 7))
  }

  pdf(paste0("Results/boxplots ",group_name,".pdf"),
      width = 11.8,
      height = 7.3)
  print(gridExtra::marrangeGrob(grobs = boxplots, ncol=4, nrow = 2))
  dev.off()
  
}

ggarrange(plotlist = volcano_plots, ncol = 1)
ggsave("Results/volcano_plots.pdf",
       width = 20,
       height = 30)


# UMAP -------------------------------------------------------------------------

umap_file <- "RDS/UMAP.RDS"
if(file.exists(umap_file)){
  umap_df <- readRDS(umap_file)
} else {
  
  samples_of_interest$File_peacoQC <- file.path("PeacoQC_Data_05_05/PeacoQC_results/fcs_files",
                                                basename(samples_of_interest$File_peacoQC))
  
  aggregate_UMAP <- make_aggregate(samples_of_interest[setdiff(rownames(samples_of_interest), only_dip_panel), ],
                                   cTotal = 50000,
                                   channels = names(markers),
                                   agg_file = "Aggregate/aggregate_umap.fcs",
                                   seed = 1,
                                   recompute = FALSE)
  set.seed(1)
  umap_res <- uwot::umap(aggregate_UMAP@exprs[, names(markers)])
  
  fsom_umap <- NewData(fsom, aggregate_UMAP)
  selection <- GetMetaclusters(fsom_umap) %in% selected_mc
  
  fsom_umap_nongr <- NewData(fsom_step2, aggregate_UMAP@exprs[selection, names(markers)])
  fsom_umap_gr <- NewData(fsom_step2_gr, aggregate_UMAP@exprs[!selection, names(markers)])
  
  labels <- rep("", nrow(aggregate_UMAP))
  labels[selection] <- cluster_label[1:107][GetClusters(fsom_umap_nongr)]
  labels[!selection] <- cluster_label[108:207][GetClusters(fsom_umap_gr)]
  
  umap_df <- data.frame(x = umap_res[,1], 
                        y = umap_res[,2], 
                        FlowSOM = labels,
                        Type = samples_of_interest[setdiff(rownames(samples_of_interest), only_dip_panel), 
                                                   "Clinical_condition_at_sampling"][aggregate_UMAP@exprs[, "File"]],
                        aggregate_UMAP@exprs)
  
  saveRDS(umap_df, umap_file)
}



library(scattermore)
umap_plots <- list()
umap_plots[["FlowSOM"]] <- ggplot(umap_df) +
  geom_scattermore(aes(x = x, y = y, col = FlowSOM), size = 2) +
  scale_color_manual(values = manual_colors) +
  ggtitle("Annotation of FlowSOM clusters") +
  theme_minimal()

umap_plots[["Type"]] <- ggplot(umap_df) +
  geom_scattermore(aes(x = x, y = y, col = Type), size = 2) +
  ggtitle("Clinical condition at sampling") +
  theme_minimal()

for(channel in names(markers)){
  umap_plots[[channel]] <- ggplot(umap_df) +
    geom_scattermore(aes_string(x = "x", y = "y", col = channel)) +
    ggtitle(markers[channel]) +
    scale_color_distiller(palette = "RdBu", direction = -1) +
    theme_minimal() +
    theme(legend.position = "none", 
          title = element_text(size = 6),
          axis.text = element_blank())
}

pdf("Results/UMAP_allcells.pdf", 
    width = 15,
    height = 6)
gridExtra::marrangeGrob(grobs = umap_plots[1:2], nrow = 1, ncol=2)
gridExtra::marrangeGrob(grobs = umap_plots[3:length(umap_plots)], ncol=9, nrow = 4)
dev.off()

# ------------------------------------------------------------------------------

lda_model <- MASS::lda(x = all_features[, 1:207], 
                      grouping = samples_of_interest$Clinical_condition_at_sampling)
lda_res <- predict(lda_model, all_features[,1:207])
ggplot(data.frame(x = lda_res$x[,1],
                  y = lda_res$x[,2],
                  class = samples_of_interest$Clinical_condition_at_sampling,
                  id = gsub(" .*", "", rownames(samples_of_interest)))) +
  geom_text(aes(x = x, y = y, col = class, label = id)) +
  theme_minimal()
ggsave("Results/lda.pdf",
       width = 11.7, height = 8.3)

