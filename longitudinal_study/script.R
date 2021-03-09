library(ggplot2)
library(flowCore)
library(FlowSOM)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(ComplexHeatmap)

source("helperFunctions.R")

if(!dir.exists("Results")) dir.create("Results")
if(!dir.exists("Results/Figures")) dir.create("Results/Figures")

## Load previous models --------------------------------------------------------

marker_info <- load_marker_excel()
DIP_channels  <- marker_info$Channel[grep("DIP", marker_info$`CyTOF panel`)]

fsom_step1 <- readRDS("model/FlowSOM_step1.RDS")

fsom_step2_gran <- readRDS("model/FlowSOM_step2_gran.RDS")
fsom_step2_gran_sub <- fsom_step2_gran
fsom_step2_gran_sub$FlowSOM$map$codes <- 
  fsom_step2_gran_sub$FlowSOM$map$codes[, DIP_channels]

fsom_step2_lymph <- readRDS("model/FlowSOM_step2_nongran.RDS")
fsom_step2_lymph <- list(FlowSOM = fsom_step2_lymph,
                         metaclustering = fsom_step2_lymph$metaclustering)
fsom_step2_lymph_sub <- fsom_step2_lymph
fsom_step2_lymph_sub$FlowSOM$map$codes <- 
  fsom_step2_lymph_sub$FlowSOM$map$codes[, DIP_channels]

## Load sample information -----------------------------------------------------

samples <- readxl::read_xlsx(path = "Metadata/Timeline_05082020_FDS.xlsx",
                             sheet = "Final_table") %>% 
  as.data.frame()
samples <- samples[samples$SampleID != "NA", ]
rownames(samples) <- samples$SampleID

only_dip_panel <- c("COVID_ICU_013_A",
                    "COVID_ICU_014_A",
                    "COVID_ICU_015_A",
                    "COVID_ICU_016_A")

samples$PatientID <- factor(samples$PatientID,
                            levels = unique(samples$PatientID[order(as.numeric(samples$overall_stay_ICU),
                                                                    samples$PatientID)]))

## Add PeacoQC information -----------------------------------------------------

PeacoQC_dir <- "PeacoQC_results_22_06"
PeacoQC_res <- read.table(file.path(PeacoQC_dir, "PeacoQC_report.txt"),
                          sep = "\t", header = TRUE, check.names = FALSE)
PeacoQC_res <- PeacoQC_res[!duplicated(PeacoQC_res$Filename, fromLast = TRUE), ]
rownames(PeacoQC_res) <- PeacoQC_res$Filename


samples$File <- file.path(PeacoQC_dir, "fcs_files",
                          paste0(sub(".fcs", "", samples$cleaned_File_Names), "_QC.fcs"))

samples$Count <- PeacoQC_res[samples$cleaned_File_Names, "Nr. Measurements"]

## Map on the FlowSOM models ---------------------------------------------------

if(file.exists("RDS/counts.RDS")){
  counts <- readRDS("RDS/counts.RDS")
  mfis <- readRDS("RDS/mfis.RDS")
} else { 
  counts <- matrix(0,
                   nrow = length(samples$SampleID),
                   ncol = fsom_step2_lymph$FlowSOM$map$nNodes + 
                          fsom_step2_gran$FlowSOM$map$nNodes,
                   dimnames = list(samples$SampleID,
                                   c(paste0("Cl", seq_len(fsom_step2_lymph$FlowSOM$map$nNodes)),
                                     paste0("grCl", seq_len(fsom_step2_gran$FlowSOM$map$nNodes)))))
  
  markers <- fsom_step2_lymph$FlowSOM$prettyColnames[fsom_step2_lymph$FlowSOM$map$colsUsed]
  mfi_names <- expand.grid(c(paste0("Cl", seq_len(fsom_step2_lymph$FlowSOM$map$nNodes)),
                             paste0("grCl", seq_len(fsom_step2_gran$FlowSOM$map$nNodes))),
                           markers)
  mfi_names <- paste0(mfi_names[,1], " ", mfi_names[,2])
  mfis <- matrix(NA,
                 nrow = length(samples$SampleID),
                 ncol = length(mfi_names),
                 dimnames = list(samples$SampleID,
                                 mfi_names))
  
  plots_lymph <- list()
  plots_gran <- list()
  sample  <- "COVID_ICU_050_D"
  for (sample in rownames(samples)) {
    message("Mapping ", sample)
    ff <- read.FCS(samples[sample, "File"])
    fsom_step1_tmp <- NewData(fsom_step1, ff)
    
    metacluster_mfis_tmp <- MetaclusterMFIs(fsom_step1_tmp)
    metacluster_mfis_tmp <- metacluster_mfis_tmp[, GetChannels(fsom_step1_tmp, 
                                                               c("CD66b", "CD45$"), 
                                                               exact = FALSE)]
    colnames(metacluster_mfis_tmp) <- GetMarkers(fsom_step1$FlowSOM, 
                                                 colnames(metacluster_mfis_tmp))
    rownames(metacluster_mfis_tmp) <-  paste0("MC",
                                              rownames(metacluster_mfis_tmp))
    
    if(sample != "COVID_ICU_040_D"){
    selected_mc_tmp <- which(metacluster_mfis_tmp[,"CD66b"] < 3 & 
                               metacluster_mfis_tmp[,"CD45"] > 3 )
    } else {
      selected_mc_tmp <- which(metacluster_mfis_tmp[,"CD66b"] < 3 & 
                                 metacluster_mfis_tmp[,"CD45"] > 4 )
    }
    selection <- GetMetaclusters(fsom_step1_tmp) %in% selected_mc_tmp
    
    
    png(paste0("Results/splitGranulocytes2/",sample,".png"),
        width = 1500, height = 1000)
    plot(ff@exprs[1:min(20000, nrow(ff)), 
                  GetChannels(ff, c(".*CD66b", ".*CD45$"), exact = FALSE)],
         pch = ".",
         col = ifelse(GetMetaclusters(fsom_step1_tmp)[1:20000] %in% selected_mc_tmp,
                      "#7bccc4",
                      "black"),
         cex = 10, cex.axis = 2, cex.lab = 2, cex.main = 2,
         xlab = "CD66b", ylab = "CD45",
         xlim = c(0,8), ylim = c(0,8),
         main = sample)
    points(metacluster_mfis_tmp[,c("CD66b", "CD45")],
           pch = 19,
           col = ifelse(seq_len(NMetaclusters(fsom_step1)) %in% selected_mc_tmp, 
                        "#2b8cbe", "#fb6a4a"),
           cex =  3)
    abline(v = 3, col = "grey", lwd = 3)
    abline(h = 3, col = "grey", lwd = 3)
    dev.off()
    
    ## _Lymphocytes --------------------------------------------------------------
    
    if (sample %in% only_dip_panel){
      fsom_lymph_tmp <- NewData(fsom_step2_lymph_sub, ff[selection,])
      fsom_lymph_tmp$prettyColnames <- fsom_step2_lymph_sub$FlowSOM$prettyColnames
      plots_lymph[[sample]] <- PlotStars(fsom_lymph_tmp,
                                           markers = DIP_channels,
                                           title = sample)
    } else {
      fsom_lymph_tmp <- NewData(fsom_step2_lymph, ff[selection,])
      fsom_lymph_tmp$prettyColnames <- fsom_step2_lymph$FlowSOM$prettyColnames
      plots_lymph[[sample]] <- PlotStars(fsom_lymph_tmp,
                                           title = sample)
    }
    
    counts_tmp <- table(GetClusters(fsom_lymph_tmp))
    counts[sample, paste0("Cl", names(counts_tmp))] <- counts_tmp
    
    if (sample %in% only_dip_panel){
      mfis_tmp <- data.frame(Cluster = paste0("Cl", seq_len(fsom_step2_lymph$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_lymph_tmp)[,DIP_channels])
      colnames(mfis_tmp) <- c("Cluster", markers[DIP_channels])
    } else {
      mfis_tmp <- data.frame(Cluster = paste0("Cl", seq_len(fsom_step2_lymph$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_lymph_tmp)[,names(markers)])
      colnames(mfis_tmp) <- c("Cluster", markers)
    }
    mfis_tmp <- tidyr::gather(mfis_tmp, "Marker", "MFI", -Cluster)
    mfis[sample, paste0(mfis_tmp$Cluster," ",mfis_tmp$Marker)] <- mfis_tmp$MFI
    
    ## _Granulocytes -------------------------------------------------------------
    
    if (sample %in% only_dip_panel){
      fsom_gran_tmp <- NewData(fsom_step2_gran_sub, ff[!selection,])
      fsom_gran_tmp$prettyColnames <- fsom_step2_gran_sub$FlowSOM$prettyColnames
      plots_gran[[sample]] <- PlotStars(fsom_gran_tmp,
                                          markers = DIP_channels,
                                          title = sample)
    } else {
      fsom_gran_tmp <- NewData(fsom_step2_gran, ff[!selection,])
      fsom_gran_tmp$prettyColnames <- fsom_step2_gran$FlowSOM$prettyColnames
      plots_gran[[sample]] <- PlotStars(fsom_gran_tmp,
                                          title = sample)
    }
    
    counts_tmp <- table(GetClusters(fsom_gran_tmp))
    counts[sample, paste0("grCl", names(counts_tmp))] <- counts_tmp
    
    if (sample %in% only_dip_panel){
      mfis_tmp <- data.frame(Cluster = paste0("grCl", seq_len(fsom_step2_gran$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_gran_tmp)[,DIP_channels])
      colnames(mfis_tmp) <- c("Cluster", markers[DIP_channels])
    } else {
      mfis_tmp <- data.frame(Cluster = paste0("grCl", seq_len(fsom_step2_gran$FlowSOM$map$nNodes)), 
                             FlowSOM::GetMFIs(fsom_gran_tmp)[,names(markers)])
      colnames(mfis_tmp) <- c("Cluster", markers)
    }
    mfis_tmp <- tidyr::gather(mfis_tmp, "Marker", "MFI", -Cluster)
    mfis[sample, paste0(mfis_tmp$Cluster," ",mfis_tmp$Marker)] <- mfis_tmp$MFI
    
  }
  
  
  saveRDS(counts, "RDS/counts.RDS")
  saveRDS(mfis, "RDS/mfis.RDS")
  
  pdf("Results/mapping_lymph.pdf",
      width = 11.7,
      height = 8.3)
  invisible(lapply(plots_lymph, plot))
  dev.off()
  
  pdf("Results/mapping_gran.pdf",
      width = 11.7,
      height = 8.3)
  invisible(lapply(plots_gran, plot))
  dev.off()

}

# Compare FlowSOM split to pathsetter ------------------------------------------

to_plot_split <- data.frame(PatientID = samples$PatientID,
                            SampleID = samples$SampleID,
                            nLymph_FlowSOM = rowSums(counts[rownames(samples),1:107]),
                            nGran_FlowSOM = rowSums(counts[rownames(samples),108:207]),
                            Lymph_PathSetter = as.numeric(samples$Lymphocytes),
                            Gran_PathSetter = as.numeric(samples$Gran),
                            Lymph_clinic = as.numeric(samples$Lymp_clin),
                            Neutr_clinic = as.numeric(samples$Neut_clin))
to_plot_split$Lymph_FlowSOM = to_plot_split$nLymph_FlowSOM / rowSums(counts[rownames(samples),]) * 100
to_plot_split$Gran_FlowSOM = to_plot_split$nGran_FlowSOM / rowSums(counts[rownames(samples),]) * 100

ggplot(to_plot_split, aes(x = Lymph_PathSetter, y = Lymph_FlowSOM)) +
  geom_point() +
  geom_text(aes(label = ifelse(abs(Lymph_FlowSOM-Lymph_PathSetter) > 20, 
                               SampleID, "")),
            nudge_y = 2,
            size = 2) +
  geom_abline(aes(intercept = 0, slope =1)) + 
  theme_minimal()

ggplot(to_plot_split, aes(x = Gran_PathSetter, y = Gran_FlowSOM)) +
  geom_point() +
  geom_text(aes(label = ifelse(abs(Gran_FlowSOM-Gran_PathSetter) > 10, 
                               SampleID, "")),
            nudge_y = 2, size = 2) +
  geom_abline(aes(intercept = 0, slope =1)) + 
  theme_minimal()

# Percentages ------------------------------------------------------------------

pctgs <- counts
pctgs <- t(apply(counts, 1, function(x) x / sum(x) * 100))


# Cluster annotation ---------------------------------------------------------------

cluster_annotation <- xlsx::read.xlsx("Metadata/FlowSOM_cluster_annotations_v4.xlsx",
                                      sheetName = "Final version",
                                      check.names = FALSE)
rownames(cluster_annotation) <- cluster_annotation$cluster

cluster_label <- cluster_annotation[c(paste0("Cl", 1:108), paste0("grCl", 1:100)), "higher level"]
cluster_label_2 <- cluster_annotation[c(paste0("Cl", 1:108), paste0("grCl", 1:100)), "second_level"]
cluster_label_final <- cluster_annotation[c(paste0("Cl", 1:108), paste0("grCl", 1:100)), "Final_annotation"]

pctgs_highlevel <- t(apply(pctgs, 1, 
                           function(sample_pctgs) tapply(sample_pctgs, 
                                                         cluster_label, 
                                                         sum)))
colnames(pctgs_highlevel) <- paste0("FlowSOM_", colnames(pctgs_highlevel))

pctgs_highlevel_2 <- t(apply(pctgs, 1, 
                             function(sample_pctgs) tapply(sample_pctgs, 
                                                           cluster_label_2, 
                                                           sum)))
colnames(pctgs_highlevel_2) <- paste0("FlowSOM_", colnames(pctgs_highlevel_2))

# Leave out duplicated columns (including CD8 NKT / CD8_NKT) and the neutrophil activation split
cols_level2_to_use <- setdiff(colnames(pctgs_highlevel_2), colnames(pctgs_highlevel))[-c(15,20,21)]

pctgs_highlevel <- cbind(pctgs_highlevel[, -c(which(colnames(pctgs_highlevel) == "FlowSOM_undefined"))],
                         "FlowSOM_lymphocytes" = rowSums(pctgs_highlevel[,c("FlowSOM_B cells",
                                                                            "FlowSOM_CD4 T cells",
                                                                            "FlowSOM_CD8 NKT",
                                                                            "FlowSOM_CD8 T cells",
                                                                            "FlowSOM_gdTcells",
                                                                            "FlowSOM_monocytes and mDC",
                                                                            "FlowSOM_NK cells",
                                                                            "FlowSOM_pDC",
                                                                            "FlowSOM_T cells DN")]),
                         pctgs_highlevel_2[, cols_level2_to_use])


comparisons <- list("Lymphocytes" =c("FlowSOM_lymphocytes", "Lymphocytes"),
                    "B cells" = c("FlowSOM_B cells", "B_Cells"),
                    "Basophils" = c("FlowSOM_Basophils", "Baso"),
                    "CD4 T" = c("FlowSOM_CD4 T cells", "CD4_T_cells"),
                    "CD8 T" = c("FlowSOM_CD8 T cells" , "CD8_T_cells"),
                    "CD8 NKT" = c("FlowSOM_CD8 NKT", "MAIT/NKT"),
                    "gd T" = c("FlowSOM_gdTcells", "GD_T"),
                    "NK"  = c("FlowSOM_NK cells","NK_Cells"),
                    "pDC" = c("FlowSOM_pDC", "pDC"),
                    "mDC" = c("FlowSOM_mDC", "mDC"),
                    "Monocytes" = c("FlowSOM_classical_monocytes", "Mono"),
                    "Neutrophils" = c("FlowSOM_neutrophil", "Neut"),
                    "Eosinophils" = c("FlowSOM_eosinophils", "Eos"),
                    "CD4_naive" = c("FlowSOM_CD4_naive", "CD4_Naive"),
                    "CD8_naive" = c("FlowSOM_CD8_naive", "CD8_Naive"))

plots <- list()

for(comp in names(comparisons)){
  x = as.numeric(samples[rownames(pctgs_highlevel), comparisons[[comp]][2]])
  y = pctgs_highlevel[,comparisons[[comp]][1]]
  plots[[comp]] <- ggplot(data.frame(SampleID = rownames(samples),
                                     x = x, 
                                     y = y), 
                          aes_string(x = x,
                                     y = y)) +
    geom_point() +
    #geom_text(aes(label = ifelse(abs(log(x/y)) > 1 & x > 0.1,
    #                             SampleID, "")),
    #          #nudge_y = 2,
    #          size = 2) +
    geom_abline(aes(intercept = 0, slope =1)) + 
    ggtitle(comp) +
    theme_minimal() +
    xlab("Pathsetter") +
    ylab("FlowSOM")
}

ggarrange(plotlist = plots)

celltypes_order <- c("FlowSOM_neutrophil",
                     "FlowSOM_eosinophils",
                     "FlowSOM_Basophils",
                     "FlowSOM_B cells",
                     "FlowSOM_gdTcells",
                     "FlowSOM_T cells DN",
                     "FlowSOM_CD4_naive",
                     "FlowSOM_CD4_effector",
                     "FlowSOM_CD4_Th1",
                     "FlowSOM_CD4_Th2",
                     "FlowSOM_CD4_Th17",
                     "FlowSOM_CD4_Tfh",
                     "FlowSOM_CD4_Treg",
                     "FlowSOM_CD8_naive",
                     "FlowSOM_CD8_activated",
                     "FlowSOM_CD8_memory",
                     "FlowSOM_CD8 NKT",
                     "FlowSOM_NK cells",
                     "FlowSOM_classical_monocytes",
                     "FlowSOM_classical_monocytes CD163+",
                     "FlowSOM_intermediate_monocytes",
                     "FlowSOM_non-classical_monocytes",
                     "FlowSOM_mDC",
                     "FlowSOM_pDC",
                     "FlowSOM_lymphocytes")

pctgs_highlevel <- pctgs_highlevel[, celltypes_order]

# tSNE plot --------------------------------------------------------------------

if(file.exists("RDS/tsne.RDS")){
  tsne <- readRDS("RDS/tsne.RDS")
  subset_tsne <- readRDS("RDS/subset_tsne.RDS")
} else { 
  set.seed(1)
  subset_tsne <- sample(seq_len(nrow(fsom_step2_lymph$FlowSOM$data)), 25000)
  tsne <- Rtsne::Rtsne(rbind(fsom_step2_gran$FlowSOM$data[subset_tsne, DIP_channels],
                             fsom_step2_lymph$FlowSOM$data[subset_tsne, DIP_channels]))
  saveRDS(subset_tsne, "RDS/subset_tsne.RDS")
  saveRDS(tsne, "RDS/tsne.RDS")
}
clusters_tsne <- c(paste0("grCl", GetClusters(fsom_step2_gran$FlowSOM)[subset_tsne]),
                   paste0("Cl", GetClusters(fsom_step2_lymph$FlowSOM)[subset_tsne]))
names(cluster_label) <- c(paste0("Cl", 1:108), paste0("grCl", 1:100))
names(cluster_label_2) <- c(paste0("Cl", 1:108), paste0("grCl", 1:100))
level1_tsne <- cluster_label[clusters_tsne]
level2_tsne <- cluster_label_2[clusters_tsne]
label_tsne <- rep("Unlabeled", length(level1_tsne))
for(level1_label in unique(level1_tsne)){
  if(paste0("FlowSOM_", level1_label) %in% celltypes_order){
    label_tsne[level1_tsne == level1_label] <- level1_label
  }
}
for(level2_label in unique(level2_tsne)){
  if(paste0("FlowSOM_", level2_label) %in% celltypes_order){
    label_tsne[level2_tsne == level2_label] <- level2_label
  }
}
label_tsne <- factor(label_tsne,
                     levels = c("Unlabeled", gsub("FlowSOM_", "", celltypes_order))) %>% 
  droplevels()

median_x <- tapply(tsne$Y[,1], label_tsne, median)
median_y <- tapply(tsne$Y[,2], label_tsne, median)

manual_colors <- c(
  "Unlabeled" = 'grey',
  "neutrophil" = "#d9c40b",
  "eosinophils" = '#ffcc99',
  "Basophils" = '#993404',
  "B cells" = '#a6bddb',
  "gdTcells" = '#ff14ef', 
  "T cells DN" = "#c288be",
  "CD4_naive" = '#d1007b',
  "CD4_effector" = '#e30237',
  "CD4_Th1" = '#ffb0de',
  "CD4_Th2" = '#dd3497',
  "CD4_Th17" = '#9c004b',
  "CD4_Tfh" = '#d15e00',
  "CD4_Treg" = '#c849d1',
  "CD8_naive" ='#88419d',
  "CD8_activated" ='#bd0af2',
  "CD8_memory" ='#c79bd4', 
  "CD8 NKT" = "#1e00ff",
  "NK cells" = '#000080',
  "classical_monocytes" = '#2ca25f',
  "classical_monocytes CD163+" = "#13ffeb",
  "intermediate_monocytes" = "#ccff00",
  "non-classical_monocytes" = "#63eb79",
  "mDC" = "#7adb35",
  "pDC" = '#eb9800')

ggplot(data.frame(tsne1 = tsne$Y[,1],
                  tsne2 = tsne$Y[,2],
                  label = label_tsne)) +
  scattermore::geom_scattermore(aes(x = tsne1, y = tsne2, color = label)) + 
  ggrepel::geom_label_repel(aes(x = x, y = y, label = label, color = label),
                           data = data.frame(x = median_x,
                                             y = median_y,
                                             label = names(median_x)),
                           segment.color = "grey", force = 20, 
                           segment.size = 0.2, point.padding = 0.5) +
  theme_minimal() +
  scale_color_manual(values = manual_colors, guide = FALSE)
ggsave(paste0("Results/Figures/210112_SuppFig6A_tSNE_populations.pdf"),
       width = 15, height = 10)

data_to_plot <- data.frame(tsne1 = tsne$Y[,1],
                           tsne2 = tsne$Y[,2],
                           rbind(fsom_step2_gran$FlowSOM$data[subset_tsne, DIP_channels],
                                 fsom_step2_lymph$FlowSOM$data[subset_tsne, DIP_channels]))
plots <- list()
for(channel in DIP_channels){
  marker <- GetMarkers(fsom_step1$FlowSOM, channel)
  plots[[marker]] <- ggplot(data_to_plot) +
    scattermore::geom_scattermore(aes_string(x = "tsne1", y = "tsne2", color = channel), pointsize = 2.2) +
    ggtitle(marker) + 
    theme_minimal() +
    scale_color_distiller(palette = "RdYlBu", direction = -1) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(color="Expression")
}
ggarrange(plotlist = plots, common.legend = TRUE)
ggsave("Results/Figures/210112_SuppFig6B_tSNE_markers.pdf",
       width = 15, height = 8)


Plot2DScatters(fsom_step2_lymph$FlowSOM,
               channelpairs = list(GetChannels(fsom_step2_lymph$FlowSOM, c("CD14","CD16$"), exact= FALSE)),
               clusters = list(which(cluster_label_2 == "classical_monocytes"),
                               which(cluster_label_2 == "intermediate_monocytes"),
                               which(cluster_label_2 == "non-classical_monocytes")),
               plotFile = paste0("Results/Figures/",date,"_monocytes.png"))

channelpair <- GetChannels(fsom_step2_lymph$FlowSOM, c("CD14","CD16$"), exact= FALSE)
data_to_plot <- data.frame(fsom_step2_lymph$FlowSOM$data[,channelpair])
p <- ggplot2::ggplot(data = data_to_plot[sample(nrow(data_to_plot), 10000), ], 
                     ggplot2::aes(x = .data$Er168Di, 
                                  y = .data$Nd148Di)) +
  ggplot2::geom_point(colour = "grey", 
                      size = 0.5) + 
  ggplot2::geom_point(data = data_to_plot[which(GetClusters(fsom_step2_lymph$FlowSOM) %in% 
                                            which(cluster_label_2 == "classical_monocytes"))[1:10000],],
                      colour = manual_colors[["classical_monocytes"]], 
                      size = 0.5) + 
  ggplot2::geom_point(data = data_to_plot[which(GetClusters(fsom_step2_lymph$FlowSOM) %in% 
                                                  which(cluster_label_2 == "intermediate_monocytes"))[1:10000],],
                      colour = manual_colors[["intermediate_monocytes"]], 
                      size = 0.5) +  
  ggplot2::geom_point(data = data_to_plot[which(GetClusters(fsom_step2_lymph$FlowSOM) %in% 
                                                  which(cluster_label_2 == "non-classical_monocytes"))[1:10000],],
                      colour = manual_colors[["non-classical_monocytes"]], 
                      size = 0.5) + 
  ggplot2::theme_classic() +
  ggplot2::xlab(GetMarkers(fsom_step2_lymph$FlowSOM, channelpair[1])) +
  ggplot2::ylab(GetMarkers(fsom_step2_lymph$FlowSOM, channelpair[2])) + 
  ggplot2::theme(legend.position = "none") +
  ggtitle("Monocyte subsets")

ggsave(plot = p,
       filename = paste0("Results/Figures/210111_monocytes_merged.png"))
# Grouping of the samples ------------------------------------------------------

neutr_lymph_ratio <- pctgs_highlevel[,"FlowSOM_neutrophil"] / 
                       pctgs_highlevel[,"FlowSOM_lymphocytes"]
names(neutr_lymph_ratio) <- rownames(pctgs_highlevel)

clus_ratio <- cut(rank(neutr_lymph_ratio), breaks = 4)
levels(clus_ratio) <- paste0("R",1:4)
clus_ratio <- factor(as.character(clus_ratio), levels = paste0("R",4:1))

# Choose clustering method -----------------------------------------------------

clusters <- clus_ratio
names(clusters) <- rownames(pctgs_highlevel)
clustername <- "Ratio"

colors <- c("R4" = "#2f4c58", 
            "R3" = "#309662",
            "R2" = "#6395ef", 
            "R1" = "#80dae8",
            "HC" = "#000000")


# Prepare samples data frame ---------------------------------------------------


samples$rank <- factor(samples$rank,
                       levels = unique(samples$rank)[
                         order(as.numeric(gsub("_.*", "", unique(samples$rank))))])
levels(samples$rank) <- c("Admission ward", "Admission ICU", "Intermediate ICU", "Discharge ICU", "Healthy control")

samples$`resp_supp_1.5` <-  as.numeric(samples[, "resp_supp_1.5"])
samples$`resp_supp_1.5`[is.na(samples$`resp_supp_1.5`)] <- 0
resp_support_levels <- c("None", "Oxygen", "HFNC", "NIV", "IMV", "IMV + Prone/ECMO/iNO")
samples$`resp_supp_1.5` <- factor(resp_support_levels[samples$`resp_supp_1.5`+1], 
                                  levels = resp_support_levels)


samples$SOFA_first_24h <- as.numeric(samples$SOFA_first_24h)
samples$Charlson_Comorbidity_Index <- as.numeric(samples$Charlson_Comorbidity_Index)
samples$days_after_discharge_ICU_overall <- as.numeric(samples$days_after_discharge_ICU_overall)
samples$age <- as.numeric(samples$age)
samples$BMI <- as.numeric(samples$BMI)


samples$BloodProfile <- clusters[rownames(samples)]
samples$Neutr_lymph_ratio <- neutr_lymph_ratio[rownames(samples)]
samples$ColorType <- as.character(samples$BloodProfile)
samples$ColorType[samples$rank == "Healthy control"] <- "HC"
samples$ColorType <- factor(samples$ColorType, levels = c(paste0("R",1:4), "HC"))

## FIGURES ---------------------------------------------------------------------

date <- "210112"

subsets <- list("all" = list(Description = "",
                             IDs = rownames(pctgs_highlevel)),
                "admission" = list(Description = "Admission samples only",
                                   IDs = samples[samples$rank == "Admission ICU", "SampleID"]),
                "discharge" = list(Description = "Discharge samples only",
                                   IDs = samples[samples$rank == "Discharge ICU", "SampleID"]),
                "highest" = list(Description = "Highest ratio per patient",
                                IDs = sapply(unique(samples$Patient), function(x){
                                  patient_samples <-  rownames(samples)[samples$Patient == x]
                                  patient_samples[which.min(clusters[patient_samples])]
                                })))

# Prepare statistics -----------------------------------------------------------

samples$`IL-8` <- samples$`IL-8_low_sens`
cytokines_APR <- c("IFN-gamma",
                   "TNF-alfa",
                   "IL-8",
                   "IL-6",
                   "IL-4",
                   "IL-2","IL-1Beta","IL-13","IL-12p70","IL-10",
                   "Eotaxin","Eotaxin_3","IP-10",
                   "MCP-1","MCP-4","MDC","MIP1a","MIP1b","TARC","GM-CSF",
                   "IL-12p40","IL-15","IL-16","IL-17","IL-1_alfa","IL-5",
                   "IL-7","TNFb","VEGF","IL1-RA","IL-18",
                   "CRP","Ddimers","Ferritin")

different_HC <- c("TNF-alfa", "IL-18", "IL-1_alfa", "MCP-1", "MIP1a")

cytokines_APR <- sort(cytokines_APR)
cytokines_APR[grep("IL", cytokines_APR)] <- cytokines_APR[grep("IL", cytokines_APR)][
  order(as.numeric(gsub("IL-([0-9]*).*", "\\1", cytokines_APR[grep("IL", cytokines_APR)])))]

to_evaluate <- list()
for(var in cytokines_APR){
  samples[[paste0("log10_", var)]] <- log10(as.numeric(samples[[var]]))
  to_evaluate[[length(to_evaluate)+1]] <-
    list(name = var,
         title = paste0(var, ifelse(var %in% different_HC, " §", "")),
         variable = paste0("log10_", var),
         subset = "all",
         test = "lmm_randomEffect")
}

clean_population_name <- function(x){
  gsub("Basophils", "basophils",
       gsub("neutrophil", "neutrophils",
            gsub("_", " ", 
                 gsub("FlowSOM_", "",x))))
}
for(population in colnames(pctgs_highlevel)){
  samples[[population]] <- pctgs_highlevel[rownames(samples), population]
  to_evaluate[[length(to_evaluate)+1]] <-
    list(name = population,
         title = clean_population_name(population),
         variable = population,
         subset = "all",
         test = "lmm_randomEffect")
}

to_evaluate_cl <- list(list(name = "Timepoint",
                            title = "Timepoint of sample",
                            variable = "rank",
                            subset = "all",
                            test = "lmm_randomEffect"),
                       list(name = "RespSupp",
                            title = "Respiratory Support",
                            variable = "resp_supp_1.5",
                            subset = "all",
                            test = "lmm_randomEffect"),
                       list(name = "SOFA",
                            title = "SOFA score first 24h",
                            variable = "SOFA_first_24h",
                            subset = "admission",
                            test = "lmm"),
                       list(name = "Revalidation",
                            title = "Days revalidation after discharge",
                            variable = "days_after_discharge_ICU_overall",
                            subset = "discharge",
                            test = "lmm"))

statistics <- data.frame(variable = character(),
                         profile1 = numeric(),
                         profile2 = numeric(),
                         pvalue = numeric(),
                         max = numeric(),
                         range = numeric(),
                         method = character())
for(evaluation in c(to_evaluate, to_evaluate_cl)){
  subset <- subsets[[evaluation$subset]]$IDs
  statistics <- statistical_tests[[evaluation$test]](
    statistics, 
    name = evaluation$name,
    variable = as.numeric(samples[subset, evaluation$variable]),
    clusters = samples[subset, "BloodProfile"],
    patientID = samples[subset, "PatientID"])
}

resp_stay <- cor.test(as.numeric(samples$max_resp_supp), 
                      as.numeric(samples$overall_stay_ICU),
                      method = "spearman")
statistics <- rbind(statistics, NA)
statistics[nrow(statistics), "variable"] <- "respiratorysupport_lengthofstay"
statistics[nrow(statistics), "pvalue"] <-  resp_stay$p.value
statistics[nrow(statistics), "method"] <- "spearman"
sofa_cor <- cor.test(as.numeric(samples$SOFA_first_24h), 
         as.numeric(samples$BloodProfile),
         method = "spearman")
statistics <- rbind(statistics, NA)
statistics[nrow(statistics), "variable"] <- "SOFA"
statistics[nrow(statistics), "pvalue"] <-  sofa_cor$p.value
statistics[nrow(statistics), "method"] <- "spearman"
discharge_cor <- cor.test(as.numeric(samples[subsets$discharge$IDs, "days_after_discharge_ICU_overall"]),
         as.numeric(samples[subsets$discharge$IDs, "BloodProfile"]),
         method = "spearman")

statistics <- rbind(statistics, NA)
statistics[nrow(statistics), "variable"] <- "discharge"
statistics[nrow(statistics), "pvalue"] <-  discharge_cor$p.value
statistics[nrow(statistics), "method"] <- "spearman"

statistics <- finalize_statistics(statistics)
statistics_filtered <- filter_statistics(statistics)

xlsx::write.xlsx(statistics[,-c(5,6)], 
                 paste0("Results/Figures/",date,"_statistics.xlsx"))

plots <- list()
for(evaluation in c(to_evaluate, to_evaluate_cl)){
  print(evaluation$name)
  p <- make_plot(title = evaluation$title,
                 variable = evaluation$variable,
                 subset = evaluation$subset)
  p <- add_stats(p, statistics_filtered, evaluation$name)
  plots[[evaluation$name]] <- p
}

# Fig 2A - Celltype boxplots ----------------------------------------------------------------------

p_celltypes <- ggarrange(plotlist = plots[colnames(pctgs_highlevel)])


# Fig 2B - Neutr/Lymph----------------------------------------------------------------------

p_neutr_lymph <- ggplot(samples[rownames(pctgs_highlevel),][order(neutr_lymph_ratio, decreasing = TRUE),], 
                        aes(x = 1:123, 
                            col = clusters[order(neutr_lymph_ratio, decreasing = TRUE)]))+
  geom_line(aes(y = FlowSOM_lymphocytes, group = 1), lwd = 0.5) +
  geom_line(aes(y = FlowSOM_neutrophil, group = 1), lwd = 0.5) +
  geom_point(aes(x = 1:123, y = FlowSOM_neutrophil,  col = ColorType, shape = rank == "Healthy control"), size = 1.5) +#alpha = rank == "Healthy control",
  geom_point(aes(x = 1:123, y = FlowSOM_lymphocytes,  col = ColorType, shape = rank == "Healthy control"), size = 1.5) + #alpha = rank == "Healthy control",, shape = 1
  ggtitle("Neutrophils and Lymphocytes") + 
  annotate("text", x = 125, y = min(pctgs_highlevel[,"FlowSOM_neutrophil"]) , label = "Neutrophils", hjust = 0, size = 4) +
  annotate("text", x = 125, y = max(pctgs_highlevel[,"FlowSOM_lymphocytes"]), label = "Lymphocytes", hjust = 0, size = 4) +
  guides(col = FALSE, shape = FALSE) +
  xlab("Samples ordered by decreasing Neutrophil/Lymphocyte ratio") +
  ylab ("Percentage") +
  coord_cartesian(xlim = c(0,125), clip = "off") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 19)) +
  theme_minimal() +
  theme(plot.margin = margin(r = 60))


test_kmeans <- kmeans(neutr_lymph_ratio, 4)
res_kmeans <- factor(test_kmeans$cluster)
p_kmeans <- ggplot(samples[rownames(pctgs_highlevel),][order(neutr_lymph_ratio, decreasing = TRUE),], 
                        aes(x = 1:123, 
                            col = res_kmeans[order(neutr_lymph_ratio, decreasing = TRUE)]))+
  geom_line(aes(y = FlowSOM_lymphocytes, group = 1), lwd = 0.5) +
  geom_line(aes(y = FlowSOM_neutrophil, group = 1), lwd = 0.5) +
  geom_point(aes(x = 1:123, y = FlowSOM_neutrophil, shape = rank == "Healthy control"), size = 1.5) +#alpha = rank == "Healthy control",
  geom_point(aes(x = 1:123, y = FlowSOM_lymphocytes, shape = rank == "Healthy control"), size = 1.5) + #alpha = rank == "Healthy control",, shape = 1
  ggtitle("Kmeans") +
  annotate("text", x = 125, y = min(pctgs_highlevel[,"FlowSOM_neutrophil"]) , label = "Neutrophils", hjust = 0, size = 4) +
  annotate("text", x = 125, y = max(pctgs_highlevel[,"FlowSOM_lymphocytes"]), label = "Lymphocytes", hjust = 0, size = 4) +
  guides(col = FALSE, shape = FALSE) +
  xlab("Samples ordered by decreasing Neutrophil/Lymphocyte ratio") +
  ylab ("Percentage") +
  coord_cartesian(xlim = c(0,125), clip = "off") +
  #scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 19)) +
  theme_minimal() +
  theme(plot.margin = margin(r = 60))

test_hclust <- hclust(dist(neutr_lymph_ratio))
res_hclust <- factor(cutree(test_hclust, 4))
p_hclust <- ggplot(samples[rownames(pctgs_highlevel),][order(neutr_lymph_ratio, decreasing = TRUE),], 
                   aes(x = 1:123, 
                       col = res_hclust[order(neutr_lymph_ratio, decreasing = TRUE)]))+
  geom_line(aes(y = FlowSOM_lymphocytes, group = 1), lwd = 0.5) +
  geom_line(aes(y = FlowSOM_neutrophil, group = 1), lwd = 0.5) +
  geom_point(aes(x = 1:123, y = FlowSOM_neutrophil, shape = rank == "Healthy control"), size = 1.5) +#alpha = rank == "Healthy control",
  geom_point(aes(x = 1:123, y = FlowSOM_lymphocytes, shape = rank == "Healthy control"), size = 1.5) + #alpha = rank == "Healthy control",, shape = 1
  ggtitle("Hierarchical clustering") +
  annotate("text", x = 125, y = min(pctgs_highlevel[,"FlowSOM_neutrophil"]) , label = "Neutrophils", hjust = 0, size = 4) +
  annotate("text", x = 125, y = max(pctgs_highlevel[,"FlowSOM_lymphocytes"]), label = "Lymphocytes", hjust = 0, size = 4) +
  guides(col = FALSE, shape = FALSE) +
  xlab("Samples ordered by decreasing Neutrophil/Lymphocyte ratio") +
  ylab ("Percentage") +
  coord_cartesian(xlim = c(0,125), clip = "off") +
  #scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 19)) +
  theme_minimal() +
  theme(plot.margin = margin(r = 60))
ggarrange(p_neutr_lymph, p_kmeans, p_hclust, nrow = 1, ncol = 3,
          labels = "AUTO")
ggsave(paste0("Results/Figures/",date,"_SuppFig7_alternativeGroupings.pdf"),
       width = 15, height = 5)
# Fig 2C - CD4/CD8 ----------------------------------------------------------------------


cd4 <- rowSums(pctgs_highlevel[,c("FlowSOM_CD4_naive",
                                  "FlowSOM_CD4_effector",
                                  "FlowSOM_CD4_Th1",
                                  "FlowSOM_CD4_Th2",
                                  "FlowSOM_CD4_Th17",
                                  "FlowSOM_CD4_Tfh",
                                  "FlowSOM_CD4_Treg")])

cd8 <- rowSums(pctgs_highlevel[,c("FlowSOM_CD8_naive",
                                  "FlowSOM_CD8_activated",
                                  "FlowSOM_CD8_memory")])


p_CD4CD8 <- ggplot(samples[rownames(pctgs_highlevel),],
                   aes(x = clusters,
                       y = cd4/cd8,
                       col = clusters,
                       shape = rank == 'Healthy control')) +
  geom_boxplot(outlier.alpha = 0, aes(group = clusters)) +
  ggbeeswarm::geom_quasirandom(aes(col = ColorType)) +
  guides(col = FALSE, shape = FALSE) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, 
                                "FALSE" = 19)) + 
  xlab("") + ylab("") +
  ggtitle("CD4/CD8 ratio")+
  theme_minimal()


# Fig 2 F - I - Clinics---------------------------------------------------------

# Supp Fig 3 - Resp Support / Revalidation / SOFA / Charlson -------------------


# Fig 2 ------------------------------------------------------------------------

p_left <- (p_celltypes / ((p_neutr_lymph | plots[["CRP"]] | plots[["Ferritin"]]) + plot_layout(nrow = 1, widths = c(2, 1, 1)))) + plot_layout(ncol = 1, heights = c(3, 1))

p_right <-  plots[["RespSupp"]] / (plots[["SOFA"]] +  coord_trans(x = "reverse")) / plots[["Timepoint"]] / plots[["Revalidation"]] 


(p_left | p_right) + plot_layout(ncol = 2, widths = c(3, 1)) + plot_annotation(tag_levels = 'A')

ggsave(paste0("Results/Figures/",
              date, "_", clustername, "_Fig2.pdf"),
       width = 20, height = 14)

# Fig 3 A - Cell types Curve fitting -------------------------------------------

library(drc)
a <- subset(samples, as.numeric(samples$intact_events)>45000)
a <- cbind(a,
           pctgs_highlevel[rownames(a),])

cur <- data.frame(val = numeric(),
                  pval = numeric())
plots_curve <- list()
for (i in colnames(pctgs_highlevel)){
  fm <- drm(a[,i]~as.numeric(a$BloodProfile), fct = L.3())
  plots_curve[[i]] <- ggplotify::as.ggplot({
    function(){
      plot(fm, log = "", 
           col="red",
           xlim=c(1,4),#3), #
           xt = 1:4,#3, #
           xtlab = paste0("R", 4:1),#3), #
           xlab = "",
           ylab = "Percentage",
           lwd = 2)
      points(as.numeric(a$BloodProfile) + rnorm(nrow(a), sd = 0.02), 
             a[,i], 
             col = colors[as.character(a$ColorType)],
             pch = c(19,1)[1 + (a$rank == "Healthy control")])
    }
  }) + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
  
  h <- summary(fm)
  cur <- rbind(cur, 
               data.frame(val = unname(fm$coefficients[3]), pval = h$coefficients[3,4]))
}
cur
rownames(cur) <- colnames(pctgs_highlevel)
cur$pval_bh <- p.adjust(cur$pval, method = "BH")
write.table(cur, file = paste0("Results/Figures/",date,"_statistics_curvefitting_cytokines.txt"), sep = "\t")

# pdf("curve_fitting_cytokines_07082020_all.pdf")
# par(mfrow=c(2,1))

p_3a <- ggarrange(plotlist = plots_curve[c("FlowSOM_classical_monocytes",
                                           "FlowSOM_eosinophils",
                                           "FlowSOM_CD8_naive",
                                           "FlowSOM_CD4_Th2",
                                           "FlowSOM_CD4_naive",
                                           "FlowSOM_CD4_Th1")],
                  ncol = 1)

# Fig 3 B - Heatmap cell types -------------------------------------------------


celltype_medians <- t(apply(pctgs_highlevel,
                            2,
                            function(x){
                              tapply(x, clusters[rownames(pctgs_highlevel)], median)
                          }))


celltype_medians_scaled <- t(apply(celltype_medians,1,function(x) (x-min(x)) / (max(x)-min(x))))

order <- rownames(dplyr::arrange(data.frame(celltype_medians_scaled),
                                 desc(R4), desc(R3), desc(R2), desc(R1)))

celltype_medians_scaled_ordered <- celltype_medians_scaled[order, ]

colors_medians <- circlize::colorRamp2(quantile(celltype_medians_scaled, c(0.01, 0.5, 0.99), 
                                                 na.rm = TRUE), 
                                        c("#e5f5f9", "#99d8c9", "#2ca25f"))
ht_medians_celltype <- ComplexHeatmap::Heatmap(celltype_medians_scaled_ordered,
                                               name = " ",
                                               cluster_rows = FALSE,
                                               cluster_columns = FALSE,
                                               column_labels = c("R4", "R3", "R2", "R1"),
                                               column_names_side = "top",
                                               column_names_rot = 0,
                                               column_names_centered = TRUE,
                                               column_title = "Relative median abundance \n of the cell populations", 
                                               column_title_side = "bottom", 
                                               column_title_gp = gpar(fontsize = 8), 
                                               col = colors_medians,
                                               heatmap_legend_param =  list(at = c(0, 1), 
                                                                            labels = c("Row min", "Row max")))

p_celltype_heatmap <- draw(ht_medians_celltype,
                           column_title = "Relative cell population abundance",
                           column_title_gp = gpar(fontface = "bold", fontsize = 12))


statistics_cell_types <- matrix(NA,
                                ncol = 3,
                                nrow = ncol(pctgs_highlevel),
                                dimnames = list(colnames(pctgs_highlevel),
                                                c("R4 vs R3", "R3 vs R2", "R2 vs R1")))
for(population in colnames(pctgs_highlevel)){
  statistics_subset <- statistics %>% dplyr::filter(population == variable)
  rownames(statistics_subset) <- paste0(statistics_subset$profile1, " vs ", statistics_subset$profile2)
  
  statistics_cell_types[population, ] <- nchar(statistics_subset[colnames(statistics_cell_types), "stars"])
  
}


ht_statistics <-Heatmap(statistics_cell_types[order,],
            name = " ",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_names_side = "top",
            column_names_rot = 0,
            column_names_centered = TRUE,
            column_title = "Differences between the NLR groups", 
            column_title_side = "bottom", 
            column_title_gp = gpar(fontsize = 8), 
            col = c("grey", "#fee8c8", "#fdbb84", "#d7301f"),
            heatmap_legend_param =  list(at = c(0, 1, 2, 3), 
                                     labels = c("NS", "*", "**", "***")))

p <- draw(ht_statistics)
ggsave(plot = grid.grabExpr(draw(p)),
       filename = "Results/Statistics_heatmap.pdf",
       width = 6, height = 12)

# Fig 3 C - Monocytes ----------------------------------------------------------


names(cluster_label_2) <- colnames(counts)

clusters_of_interest <- c("Cl21", "Cl11", "Cl44", "Cl31", "Cl3", "Cl12", "Cl43",
                          "Cl33", "Cl34", "Cl13", "Cl63", "Cl41", "Cl24", "Cl32",
                          "Cl22", "Cl53","Cl23", "Cl1", "Cl2","Cl108")

markers_of_interest <- c("CD45$", "HLA-DR", "CD16$", "CD11c", "CD38", "CD14", "IL-7Ra", "CD163", "CD45RA")

mfis <- GetMFIs(fsom_step2_lymph$FlowSOM)[as.numeric(gsub("Cl", "", clusters_of_interest)), 
                                                GetChannels(fsom_step2_lymph$FlowSOM, markers_of_interest, exact = FALSE)]
rownames(mfis) <- clusters_of_interest
colnames(mfis) <- markers_of_interest


patient_overview <- t(apply(pctgs[,clusters_of_interest], 2, 
                            function(cluster_pctgs) 
                              tapply(cluster_pctgs, clusters, median)))


patient_overview_scaled <- t(apply(patient_overview,1,function(x) (x-min(x)) / (max(x)-min(x))))

order <- rownames(dplyr::arrange(data.frame(patient_overview_scaled == 1), 
                                 desc(R4), desc(R3), desc(R2), desc(R1))) 

mfis_ordered <- mfis[order, ]
patient_overview_scaled_ordered <- patient_overview_scaled[order, ]

colors_patients <- circlize::colorRamp2(quantile(patient_overview_scaled, c(0.01, 0.5, 0.99), 
                                                 na.rm = TRUE), 
                                        c("#e5f5f9", "#99d8c9", "#2ca25f"))
ht_patients <- Heatmap(patient_overview_scaled_ordered,
                       name = " ",
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       row_labels = rep("", nrow(patient_overview_scaled_ordered)),
                       column_names_side = "top",
                       column_title = "Relative median distribution\n(% of non-granulocytes)", 
                       column_title_side = "bottom", 
                       column_title_gp = gpar(fontsize = 8), 
                       col = colors_patients,
                       heatmap_legend_param =  list(at = c(0, 1), 
                                                    labels = c("Row min", "Row max"))#,
                       #show_heatmap_legend = FALSE
                       )


#colors_mfi <- circlize::colorRamp2(quantile(mfis, seq(0.01, 0.99, length.out = 7), 
colors_mfi <- circlize::colorRamp2(quantile(mfis, c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.99), 
                                            na.rm = TRUE), 
                                   rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))
colnames(mfis_ordered) <- gsub("\\$", "", colnames(mfis_ordered))
ht_mfi <- Heatmap(mfis_ordered,
                  cluster_columns = FALSE,
                  name = "Expression",
                  row_labels =  paste0(rownames(mfis_ordered), ": ", 
                                       gsub("_",  " ", cluster_label_2[rownames(mfis_ordered)])),
                  column_names_side = "top",
                  column_title = "Median expression level\n(arcsinh transformed \nwith cofactor 5)", 
                  column_title_side = "bottom", 
                  column_title_gp = gpar(fontsize = 8),
                  col = colors_mfi)#,
                  #show_heatmap_legend = FALSE)

p_monocytes <- draw(`+.AdditiveUnit`(ht_patients,
                                     ht_mfi),
                    column_title = "Differential monocyte clusters", 
                    column_title_gp = gpar(fontface = "bold", fontsize = 12), 
                    auto_adjust = FALSE)

ggsave(plot = grid.grabExpr(draw(p_monocytes)),
       filename = "Results/Figures/Monocytes_heatmap.pdf")

# Fig 3 ------------------------------------------------------------------------


ggarrange(p_3a, 
          grid.grabExpr(draw(p_celltype_heatmap)), 
          grid.grabExpr(draw(p_monocytes)),
          #grid.grabExpr(draw(ComplexHeatmap::Legend(at = c(0,1), 
          #                       labels = c("Row min", "Row max"), 
          #                       col_fun = colors_patients))),
          labels = "AUTO",
          widths = c(1,1,2),
          ncol = 3)


ggsave(paste0("Results/Figures/",
              date, "_", clustername, "_Fig3.pdf"),
       width = 25, height = 12)



# Fig 4B - Cytokine heatmap ----------------------------------------------------

cytokines_to_plot <- cytokines_APR[cytokines_APR %in% statistics_filtered$variable]
cytokines_to_plot <- cytokines_to_plot[-c(1,3)] # Remove CRP and ferritin

names(to_evaluate) <- lapply(to_evaluate, function(x) x$name)


means <- data.frame(variable =character(),
                    profile1 = numeric(),
                    profile2 = numeric(),
                    profile3 = numeric(),
                    profile4 = numeric())

for(evaluation in to_evaluate[cytokines_to_plot]){
  
  means <- add_lmm_randomEffect_means(
    means, 
    name = evaluation$name,
    variable = samples[, evaluation$variable],
    clusters = samples[, "BloodProfile"],
    patientID = samples[, "PatientID"])
}
rownames(means) <- means$variable
means <- means[2:5]

means_scaled <- t(apply(means,1,function(x) (x-min(x)) / (max(x)-min(x))))

order <- rownames(dplyr::arrange(data.frame(means_scaled),#data.frame(means_scaled == 1),
                                 desc(profile1), desc(profile2), desc(profile3), desc(profile4)))

means_scaled_ordered <- means_scaled[order, ]

colors_patients <- circlize::colorRamp2(quantile(means_scaled, c(0.01, 0.5, 0.99), 
                                                 na.rm = TRUE), 
                                        c("#e5f5f9", "#99d8c9", "#2ca25f"))
ht_means <- Heatmap(means_scaled_ordered,
                    name = " ",
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_labels = c("R4", "R3", "R2", "R1"),
                    column_names_side = "top",
                    column_names_rot = 0,
                    column_names_centered = TRUE,
                    column_title = "Relative mean cytokine values\n(estimated by lmm with random patient effect)", 
                    column_title_side = "bottom", 
                    column_title_gp = gpar(fontsize = 8), 
                    col = colors_patients,
                    heatmap_legend_param =  list(at = c(0, 1), 
                                                 labels = c("Row min", "Row max")))


# Fig 4A - Cytokines -----------------------------------------------------------


p_4a <- gridExtra::grid.arrange(grobs = plots[order],
                                layout_matrix = matrix(c(1:20), ncol = 5, byrow = TRUE))


# Fig 4 ------------------------------------------------------------------------


ggarrange(p_4a, 
          ggparagraph(text=" ", face = "italic", size = 16, color = "white"), 
          grid.grabExpr(draw(ht_means)),
          labels = c("A", "B", ""),
          widths = c(3,0.1,1),
          nrow = 1)


ggsave(paste0("Results/Figures/",
              date, "_", clustername, "_Fig4.pdf"),
       width = 20, height = 12)




# Supp Fig 1 - Heatmap Bloodprofile --------------------------------------------


annotation <- data.frame("Ratio group" = clusters, check.names = FALSE)
rownames(annotation) <- rownames(pctgs_highlevel)

dev.off()
pdf(paste0("Results/Figures/",
           date, "_", clustername, "_SuppFig1_Heatmap_Bloodprofiles.pdf"),
    width = 20, height = 8)
pheatmap::pheatmap(t(pctgs_highlevel[order(neutr_lymph_ratio, decreasing = TRUE), ]),
                   gaps_col = cumsum(table(clusters)),
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   breaks = (1:100)^3 / (100)^3*100,
                   annotation_col = annotation,
                   annotation_colors = list("Ratio group" = colors[levels(clusters)]),
                   fontsize_col = 6)
dev.off()

# Supp Fig 2 - Counts ----------------------------------------------------------

samples$volume <- as.numeric(samples$intact_events) / (as.numeric(samples$WBC_clin) * 10^9)

abs_counts <- data.frame("LymphocytesCyTOF" =  (((pctgs_highlevel[,"FlowSOM_lymphocytes"]/100) *
                                                   as.numeric(samples[rownames(pctgs_highlevel), "intact_events"]))/ 
                                                  samples[rownames(pctgs_highlevel), "volume"])/10^9,
                         "LymphocytesClinic" = as.numeric(samples[rownames(pctgs_highlevel), "Lymp_clin"]),
                         "NeutrophilsCyTOF" =  (((pctgs_highlevel[,"FlowSOM_neutrophil"]/100) *
                                                   as.numeric(samples[rownames(pctgs_highlevel), "intact_events"]))/ 
                                                  samples[rownames(pctgs_highlevel), "volume"])/10^9,
                         "NeutrophilsClinic" = as.numeric(samples[rownames(pctgs_highlevel), "Neut_clin"]),
                         "clusters" = clusters,
                         "rank" = samples[rownames(pctgs_highlevel), "rank"],
                         "ColorType" = samples[rownames(pctgs_highlevel), "ColorType"])
abs_counts$RatioCyTOF <- abs_counts$NeutrophilsCyTOF / abs_counts$LymphocytesCyTOF
abs_counts$RatioClinic <- abs_counts$NeutrophilsClinic / abs_counts$LymphocytesClinic

plots <- list()

y_labels <- c("LymphocytesCyTOF" = "10^9 cells / l",
              "NeutrophilsCyTOF" = "10^9 cells / l",
              "RatioCyTOF" = "Neutr. / Lymph.",
              "LymphocytesClinic" = "10^9 cells / l",
              "NeutrophilsClinic" = "10^9 cells / l",
              "RatioClinic" = "Neutr. / Lymph.")
for(to_plot in c("LymphocytesCyTOF",
                 "NeutrophilsCyTOF",
                 "RatioCyTOF",
                 "LymphocytesClinic",
                 "NeutrophilsClinic",
                 "RatioClinic")){
  plots[[length(plots) + 1]] <- 
    ggplot(abs_counts,
           aes_string(x = "clusters", 
                      y = to_plot,
                      col = "clusters",
                      shape = "rank == 'Healthy control'")) +
    geom_boxplot(outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom(aes(col=ColorType)) +
    ggtitle(gsub("C", " measured by C", to_plot)) + 
    xlab("") + ylab(y_labels[to_plot]) +
    guides(col = FALSE, shape = FALSE) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c("TRUE" = 1, 
                                  "FALSE" = 19)) +
    theme_minimal()
}

for(to_plot in c("Lymphocytes", "Neutrophils", "Ratio")){
  cor <- cor(abs_counts[,paste0(to_plot, "CyTOF")],
             abs_counts[,paste0(to_plot, "Clinic")],
             use = "complete.obs")
  
  plots [[length(plots) + 1]] <-
    ggplot(abs_counts,
           aes_string(x = paste0(to_plot, "CyTOF"),
                      y = paste0(to_plot, "Clinic"),
                      col = "clusters",
                      shape = "rank == '11_healthy_control'")) +
    geom_point() +
    geom_abline(slope = 1) +
    ggtitle(paste0(to_plot," (r = ",round(cor,3),")")) + 
    xlab("Measured by CyTOF") + ylab("Measured by Clinic") +
    guides(col = FALSE, shape = FALSE) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c("TRUE" = 1, 
                                  "FALSE" = 19)) +
    theme_minimal()
}

pdf(paste0("Results/Figures/",
           date, "_", clustername, "_SuppFig2_Boxplots_CountsPerVolume.pdf"),
    width = 14, height = 8)
gridExtra::grid.arrange(grobs = plots, layout_matrix = matrix(1:9, ncol = 3))

dev.off()

# Supp Fig 3 - patient verloop  ------------------------------------------------------------------


samples$Date_of_sample_event <- strptime(samples$Date_of_sample_event, 
                                         format = "%d/%m/%Y__%H:%M:%S")

patients <- grep("HC", unique(samples$PatientID), value = TRUE, invert = TRUE)
blood_profiles <- list()
pattern_type <- list()
for(patient in patients){
  samples_patient <- samples[samples$PatientID == patient, ]
  samples_patient <- samples_patient[order(samples_patient$Date_of_sample_event), ]
  df <- cbind(samples_patient,
              Order = seq_len(nrow(samples_patient)))
  
  diff <-   diff(as.numeric(df$BloodProfile[df$rank != "Admission ward"]))
  
  pattern_type[[patient]] <- ifelse(all(diff == 0),
                                    "No change",
                                    ifelse(all(diff <= 0),
                                           "Increasing",
                                           ifelse(all(diff >= 0),
                                                  "Improving",
                                                  "Variable")))
  df$pattern_type <-  pattern_type[[patient]]
  
  blood_profiles[[patient]] <- df
}

pattern_type <- unlist(pattern_type)
pattern_type <- factor(pattern_type,
                       levels = c("Improving", "No change", "Variable", "Increasing"))
table(pattern_type)


blood_profiles_all <- do.call(rbind, blood_profiles)
rownames(blood_profiles_all) <- gsub(".*\\.", "", rownames(blood_profiles_all))

blood_profiles_all$PatientID <- factor(as.character(blood_profiles_all$PatientID),
                                       levels = names(sort(pattern_type)))


blood_profiles_all$BloodProfile <- factor(as.character(blood_profiles_all$BloodProfile),
                                          levels = paste0("R",1:4))

p_patientChange <- ggplot(blood_profiles_all,
                          aes(x = as.numeric(Days_from_ad_ICU),
                              y = BloodProfile)) +
  geom_line(aes(group = PatientID)) + 
  geom_point(aes(col = rank)) +
  theme_minimal() +
  facet_wrap(~ PatientID) +
  xlab("Days since admission to ICU") +
  scale_color_manual(values = c("Admission ward" = "#BBBBBB",
                                "Admission ICU" = "#FF0000",
                                "Intermediate ICU" = "#FFA500",
                                "Discharge ICU" = "#006400"))


pdf(paste0("Results/Figures/",
           date, "_", clustername, "_SuppFig3_PatientChange.pdf"),
    width = 14, height = 8)
print(p_patientChange)
dev.off()


# Supp Fig 4 A - BMI etc ---------------------------------------------------------

to_evaluate <- list(list(name = "Age",
                         title = "Age",
                         variable = "age",
                         subset = "highest"),
                    list(name = "BMI_lowest",
                         title = "BMI",
                         variable = "BMI",
                         subset = "highest"),
                    list(name = "BMI_discharge",
                         title = "BMI",
                         variable = "BMI",
                         subset = "discharge"),
                    list(name = "Comorbidity",
                         title = "Charlson Comorbidity Index",
                         variable = "Charlson_Comorbidity_Index",
                         subset = "admission"))

samples$colorType <- as.character(samples$ColorType)
supp_4A <- list()
for(evaluation in to_evaluate){
  print(evaluation$name)
  p <- make_plot(title = evaluation$title,
                 variable = evaluation$variable,
                 subset = evaluation$subset)
  supp_4A[[evaluation$name]] <- p
}


# Supp Fig 4 B -----------------------------------------------------------------
data_supp_4b <- data.frame("Age"  = samples[subsets$highest$IDs, "age"],
                           "BMI" = samples[subsets$highest$IDs, "BMI"],
                           "NLR" = samples[subsets$highest$IDs, "Neutr_lymph_ratio"],
                           "Group" = samples[subsets$highest$IDs, "colorType"])
data_supp_4b <- data_supp_4b[!is.na(data_supp_4b$Age),]

cor_nlr_age <- cor(data_supp_4b$Age, data_supp_4b$NLR, use = "complete.obs")
supp4b_age <- ggplot(data_supp_4b) +
  geom_point(aes(x = NLR, y = Age, color = Group))+
  scale_color_manual(values = colors) +
  theme_minimal() +
  scale_x_reverse()


supp4b_bmi <- ggplot(data_supp_4b) +
  geom_point(aes(x = NLR, y = BMI, color = Group))+
  scale_color_manual(values = colors) +
  theme_minimal() +
  scale_x_reverse()

cor_nlr_bmi <- cor(data_supp_4b$BMI, data_supp_4b$NLR, use = "complete.obs")

# samples_unique <- 
#   t(sapply(unique(as.character(samples$PatientID)),
#            function(id) c(id,
#                           samples[samples$PatientID == id, "max_resp_supp"][1],
#                           samples[samples$PatientID == id, "BMI"][1],
#                           samples[samples$PatientID == id, "age"][1],
#                           samples[samples$PatientID == id, "gender"][1],
#                           samples[samples$PatientID == id, "SOFA_first_24h"][1],
#                           samples[samples$PatientID == id, "Charlson_Comorbidity_Index"][1]))) %>% 
#   data.frame() %>% 
#   magrittr::set_colnames(c("PatientID", "max_resp_supp", "BMI", 
#                            "age", "gender", "SOFA_first_24h", "Charlson_Comorbidity_Index"))
# 
# for(i in c( "max_resp_supp", "BMI", "SOFA_first_24h", "Charlson_Comorbidity_Index")){
#   samples_unique[[i]] <- as.numeric(samples_unique[[i]])
# }
# for(i in c( "max_resp_supp", "SOFA_first_24h", "Charlson_Comorbidity_Index")){
#   samples_unique[[i]] <- factor(samples_unique[[i]])
# }
# 
# samples_unique <- samples_unique[!is.na(samples_unique$max_resp_supp), ]
# levels(samples_unique$max_resp_supp) <- levels(samples$resp_supp_1.5)[2:6]
# 
# p_BMI_RespSupp <- ggplot(samples_unique,
#                          aes(x = factor(max_resp_supp),
#                              y = BMI)) +
#   geom_boxplot(outlier.alpha = 0) +
#   ggbeeswarm::geom_quasirandom() +
#   guides(col = FALSE) +
#   xlab("Maximal Respiratory Support") + ylab("BMI") +
#   ggtitle("BMI vs Maximal Respiratory Support")+
#   theme_minimal()
# 
# p_BMI_SOFA <- ggplot(samples_unique,
#                      aes(x = factor(SOFA_first_24h),
#                          y = BMI)) +
#   geom_boxplot(outlier.alpha = 0) +
#   ggbeeswarm::geom_quasirandom() +
#   guides(col = FALSE) +
#   xlab("SOFA first 24h") + ylab("BMI") +
#   ggtitle("BMI vs SOFA")+
#   theme_minimal()
# 
# p_Comorbidity_RespSupp <- ggplot(samples_unique,
#                                  aes(x = factor(max_resp_supp),
#                                      y = as.numeric(Charlson_Comorbidity_Index))) +
#   geom_boxplot(outlier.alpha = 0, aes(group = max_resp_supp)) +
#   ggbeeswarm::geom_quasirandom() +
#   guides(col = FALSE) +
#   xlab("Maximal Respiratory Support") + ylab("Charlson Comorbidity Index") +
#   ggtitle("Comorbidity vs Maximal Respiratory Support")+
#   theme_minimal()


# Supp Fig 4 C - Corticosteroids ------------------------------------------------------------

samples$`Anti-IL-1`[samples$`Anti-IL-1` == 1] <- "Anti-IL-1"
samples$`Anti-IL-1`[samples$`Anti-IL-1` == 0] <- ""
samples$`Anti-IL-6`[samples$`Anti-IL-6` == 1] <- "Anti-IL-6"
samples$`Anti-IL-6`[samples$`Anti-IL-6` == 0] <- ""

samples$ILtreatment <- ifelse(samples$`Anti-IL-1` == "Anti-IL-1", 
                              "Anti-IL-1",
                              ifelse(samples$`Anti-IL-6` == "Anti-IL-6",
                                     "Anti-IL-6",
                                     "No IL-1 or IL-6 treatment"))

p_steroids_discharge <- 
  ggplot(samples[subsets$discharge$IDs,],
         aes(x = ifelse(CS == "1", 
                        "Treated with corticosteroids", 
                        "No corticosteroids"),
             y = factor(as.character(BloodProfile), 
                        levels = paste0("R",1:4)))) +
  geom_boxplot(outlier.alpha = 0,
               aes(group =  CS)) +
  ggbeeswarm::geom_quasirandom(aes(col = BloodProfile,
                                   shape = ILtreatment), 
                               size = 3) +
  ggtitle("Corticosteroid treatment") + 
  xlab("Discharge samples only") + ylab("") +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("No IL-1 or IL-6 treatment" = 19,
                                "Anti-IL-1" = 10,
                                "Anti-IL-6" = 13),
                     name = "Additional treatment") +
  theme_minimal()


# Supp Fig 4 -------------------------------------------------------------------

# ((( supp_4A$BMI_lowest | supp_4A$Age) / 
#    (p_Comorbidity_RespSupp | supp_4A$Comorbidity))| p_steroids_discharge )  + 
((( supp_4A$BMI_lowest | supp_4A$Age) / 
    (supp4b_bmi | supp4b_age))| p_steroids_discharge )  + 
  plot_layout(ncol = 2, widths = c(3, 1)) +
  plot_annotation(tag_levels = 'A')


ggsave(paste0("Results/Figures/",
           date, "_", clustername, "_SuppFig4.pdf"),
    width = 20, height = 12)

# Supp Fig 5 -------------------------------------------------------------------

mfis_lymph <- GetMFIs(fsom_step2_lymph, colsUsed = TRUE, prettyColnames = TRUE)
rownames(mfis_lymph) <- paste0("Cl", rownames(mfis_lymph))
colnames(mfis_lymph) <- gsub(".*_(.*) <.*", "\\1", colnames(mfis_lymph))

mfis_gran <- GetMFIs(fsom_step2_gran, colsUsed = TRUE, prettyColnames = TRUE)
rownames(mfis_gran) <- paste0("grCl", rownames(mfis_gran))


colors_mfi <- circlize::colorRamp2(c(0:6),
#quantile(mfis_gran, seq(0.01, 0.99, length.out = 7),
#                                             na.rm = TRUE), 
                                   rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))

pdf(paste0("Results/Figures/",date,"_",clustername,"_SuppFig5_Heatmaps.pdf"),
    width = 12, height = 20)
Heatmap(mfis_gran,
        cluster_columns = FALSE,
        name = "Expression",
        row_labels =  paste0(rownames(mfis_gran), ": ", 
                             gsub("_",  " ", cluster_label_2[rownames(mfis_gran)])),
        column_names_side = "top",
        column_title = "Median expression level\n(arcsinh transformed \nwith cofactor 5)", 
        column_title_side = "bottom", 
        column_title_gp = gpar(fontsize = 8),
        col = colors_mfi)

Heatmap(mfis_lymph,
        cluster_columns = FALSE,
        name = "Expression",
        row_labels =  paste0(rownames(mfis_lymph), ": ", 
                             gsub("_",  " ", cluster_label_2[rownames(mfis_lymph)])),
        column_names_side = "top",
        column_title = "Median expression level\n(arcsinh transformed \nwith cofactor 5)", 
        column_title_side = "bottom", 
        column_title_gp = gpar(fontsize = 8),
        col = colors_mfi)
dev.off()


## Export data table

data_export <- data.frame(PatientID = rownames(pctgs_highlevel),
                          pctgs_highlevel, 
                          samples[rownames(pctgs_highlevel), cytokines_APR],
                          Group = samples[rownames(pctgs_highlevel), "BloodProfile"],
                          check.names = FALSE)

writexl::write_xlsx(data_export,
                    paste0(date,"_data_export.xlsx"))


###
ggplot(data.frame("Maximal respiratory support" = samples$max_resp_supp[samples$rank != "11_healthy_control"],
                  "Length of stay ICU" = as.numeric(samples$overall_stay_ICU[samples$rank != "11_healthy_control"]),
                  check.names = FALSE)) +
  geom_point(aes(x = `Maximal respiratory support`, 
                 y = `Length of stay ICU`)) + 
  theme_minimal()
ggsave("Results/Figures/maxRespSupp_vs_overallStayICU.pdf")

###


mfis <- readRDS("RDS/mfis.RDS")
mfis_hladr  <- mfis[,grep("HLA-DR", colnames(mfis))]
colnames(mfis_hladr) <- gsub(" .*", "", colnames(mfis_hladr))

mfis_hladr <- tidyr::pivot_longer(data.frame(Patient = rownames(mfis_hladr),
                                             mfis_hladr), 
                                  names_to = "Cluster", 
                                  values_to = "MMI",
                                  cols = colnames(mfis_hladr))
mfis_hladr$Celltype <- cluster_label_2[mfis_hladr$Cluster]
mfis_hladr$BloodProfile <- samples[mfis_hladr$Patient, "BloodProfile"]


ggplot(mfis_hladr %>% 
         dplyr::filter(Celltype %in% c("classical_monocytes", "intermediate_monocytes", "non-classical_monocytes", "mDC", "pDC"))) +
  geom_boxplot(aes(x = BloodProfile, y = MMI), outlier.size = 0) +
  ggbeeswarm::geom_quasirandom(aes(x = BloodProfile, y = MMI, color = Celltype)) +
  facet_wrap("~ Cluster", nrow = 1) +
  scale_color_manual(values = manual_colors) + 
  theme_minimal() +
  xlab("Immune type")

ggsave("Results/Figures/210112_HLADR_MMIs.pdf",
       width = 20, height = 8)

mfis_hladr <- mfis_hladr %>%
  dplyr::group_by(Celltype, Patient, BloodProfile) %>%
  dplyr::summarise(MMI = median(MMI))

ggplot(mfis_hladr %>% 
         dplyr::filter(Celltype %in% c("classical_monocytes", "intermediate_monocytes", "non-classical_monocytes", "mDC", "pDC"))) +
  geom_boxplot(aes(x = BloodProfile, y = MMI), outlier.size = 0) +
  ggbeeswarm::geom_quasirandom(aes(x = BloodProfile, y = MMI, color = Celltype)) +
  facet_wrap("~ Celltype", nrow = 1) +
  scale_color_manual(values = manual_colors) + 
  theme_minimal() +
  xlab("Immune type")


ggsave("Results/Figures/210112_HLADR_median_of_MMIs.pdf",
       width = 15, height = 8)


plot(pctgs[,c("Cl11","grCl63")], xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "lightgrey")
# tSNE test -------------------------------

data_tsne <- data.frame(Ratio = samples[rownames(pctgs_highlevel),"BloodProfile"],
                        PatientID = samples[rownames(pctgs_highlevel),"PatientID"],
                        pctgs_highlevel, 
                        apply(samples[rownames(pctgs_highlevel),cytokines_APR], 2, as.numeric),
                        check.names = FALSE)
rownames(data_tsne) <- rownames(pctgs_highlevel)
data_tsne  <- na.omit(data_tsne)
set.seed(1)
res_tsne <- Rtsne::Rtsne(scale(data_tsne[,-c(1,2)]), perplexity = 10)
data_tsne$tsne1 <- res_tsne$Y[,1]
data_tsne$tsne2 <- res_tsne$Y[,2]

tsne_plots <- list()
for(feature in colnames(data_tsne)[c(1,2)]){
  tsne_plots[[feature]] <- ggplot(data_tsne) +
    geom_point(aes_string(x = "tsne1", y = "tsne2", col = paste0("`",feature,"`"))) +
    #theme_minimal() +
    ggtitle(feature)
}
for(feature in colnames(data_tsne)[-c(1,2)]){
  tsne_plots[[feature]] <- ggplot(data_tsne) +
    geom_point(aes_string(x = "tsne1", y = "tsne2", col = paste0("`",feature,"`"))) +
    #theme_minimal() +
    scale_color_distiller(palette = "RdYlBu") +
    ggtitle(feature)
}

pdf("Results/tsne_test.pdf", width = 20, height = 10)
ggarrange(plotlist = tsne_plots, ncol = 3, nrow = 2)
dev.off()


tsne_plots <- list()
for(feature in colnames(data_tsne)[c(1,2)]){
  tsne_plots[[feature]] <- ggplot(data_tsne) +
    geom_point(aes_string(x = "FlowSOM_neutrophil", y = "FlowSOM_lymphocytes", col = paste0("`",feature,"`"))) +
    #theme_minimal() +
    ggtitle(feature)
}
for(feature in colnames(data_tsne)[-c(1,2)]){
  tsne_plots[[feature]] <- ggplot(data_tsne) +
    geom_point(aes_string(x = "FlowSOM_neutrophil", y = "FlowSOM_lymphocytes", col = paste0("`",feature,"`"))) +
    #theme_minimal() +
    scale_color_distiller(palette = "RdYlBu") +
    ggtitle(feature)
}

pdf("Results/ratio_test.pdf", width = 20, height = 10)
ggarrange(plotlist = tsne_plots, ncol = 3, nrow = 2)
dev.off()


# SOFA score -----------------

sofa_test <- data.frame(SOFA = samples$SOFA_first_24h,
                        RatioGroup = samples$BloodProfile,
                        PatientID = samples$PatientID,
                        Timepoint = samples$rank,
                        Days_from_ad_ICU = as.numeric(samples$Days_from_ad_ICU),
                        Ratio = samples$Neutr_lymph_ratio)


ggplot(sofa_test) +
  geom_point(aes(x = Days_from_ad_ICU , y = Ratio, col = RatioGroup,
                 size = Timepoint, shape = Timepoint))  +
  geom_line(aes(x = Days_from_ad_ICU, y = Ratio, group = PatientID)) +
  facet_wrap(~ SOFA) +
  scale_size_manual(values = c("Admission ICU" = 2,
                               "Admission ward" = 1,
                               "Intermediate ICU" = 1,
                               "Discharge ICU" = 1,
                               "Healthy control" = 1)) +
  ggtitle("SOFA score groups") +
  scale_color_manual(values = colors)
ggsave("Results/Sofa_exploration.pdf",
       width = 10, height = 6)

###

Plot2DScatters(fsom_step2_lymph$FlowSOM, list(c("CD11c", "CD38"),
                                              c("CD14", "CD16"),
                                              c("CD14", "IL-3R"),
                                              c("CD14", "CD45RA")),
               clusters = list(c(11,53,63,41,32,54,31,13,24,22,21,44,3,12,43,33,34,23,1,2)), 
               plotFile = "Monocyte_scatters_all.png",
               maxPoints = 50000)

# ___Old plots___---------------------------------------------------------------



# Fig 5 - Cytokines ------------------------------------------------------------

p_neutr_lymph <- ggplot(cbind(samples[rownames(pctgs_highlevel),],
                              pctgs_highlevel)[order(clusters),],
                        aes(x = 1:123, 
                            col = clusters[order(clusters)])) +
  geom_line(aes(y = FlowSOM_lymphocytes, group = 1), lwd = 1) +
  geom_line(aes(y = FlowSOM_neutrophil, group = 1), lwd = 1) +
  ggtitle("Neutrophils and Lymphocytes") + 
  annotate("text", x = 125, y = 56.63363 , label = "Neutrophils", hjust = 0, size = 4) +
  annotate("text", x = 125, y = 33.634208, label = "Lymphocytes", hjust = 0, size = 4) +
  guides(col = FALSE) +
  xlab("Samples ordered by blood profile 1-4") +
  ylab ("Percentage") +
  coord_cartesian(xlim = c(0,125), clip = "off") +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(plot.margin = margin(r = 60))



pdf(paste0("Results/Figures/",
           date, "_", clustername, "_Fig5_boxplotsCytokines.pdf"),
    width = 22, height = 12)
# Significant cytokines
cytokines_to_plot <- cytokines_APR[cytokines_APR %in% statistics_filtered$variable]

gridExtra::grid.arrange(grobs = c(list(p_neutr_lymph = p_neutr_lymph),
                                  plots[cytokines_to_plot]),
                        #layout_matrix = matrix(c(1,1, 2:5, 1, 1, 6:21), ncol = 6, byrow = TRUE))
                        layout_matrix = matrix(c(1,1, 2:23), ncol = 6, byrow = TRUE))
dev.off()

# Fig 4 - Cytokine heatmap -----------------------------------------------------


pdf(paste0("Results/Figures/",
           date, "_", clustername, "_Fig4_Heatmap_Cytokines.pdf"),
    width = 5, height = 6)
ht_means
dev.off()



# Boxplots_Celltypes -----------------------------------------------------------

pctgs_highlevel_long <- tidyr::gather(data.frame(SampleID = samples[rownames(pctgs_highlevel), "SampleID"],
                                                 Cluster = clusters,
                                                 pctgs_highlevel,
                                                 check.names = FALSE),
                                      "Celltype",
                                      "Pctgs",
                                      -SampleID, -Cluster) 

# Clean up names a bit
pctgs_highlevel_long$Celltype <- factor(gsub("Basophils", "basophils",
                                             gsub("neutrophil", "neutrophils",
                                                  gsub("_", " ", 
                                                       gsub("FlowSOM_", "", pctgs_highlevel_long$Celltype)))),
                                        levels = gsub("Basophils", "basophils",
                                                      gsub("neutrophil", "neutrophils",
                                                           gsub("_", " ", 
                                                                gsub("FlowSOM_", "", colnames(pctgs_highlevel))))))


p_celltypes <- ggplot(pctgs_highlevel_long, aes(x = Cluster,
                                                y = Pctgs,
                                                col = Cluster)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  facet_wrap(~Celltype, scales = "free") +  
  guides(col = FALSE) +
  theme_minimal() +
  xlab("")  + ylab("Percentage") +
  scale_color_manual(values = colors)

ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_Boxplots_Celltypes.pdf"),
       width = 12, height = 8)


## STATISTICS ------------------------------------------------------------------


# Boxplots_ClinicalOutcomes ----------------------------------------------------

to_evaluate <- list(list(name = "Timepoint",
                         title = "Timepoint of sample",
                         variable = "rank",
                         subset = "all",
                         test = "fisher"),
                    list(name = "SOFA",
                         title = "SOFA score first 24h",
                         variable = "SOFA_first_24h",
                         subset = "admission",
                         test = "fisher"),
                    list(name = "Comorbidity",
                         title = "Charlson Comorbidity Index",
                         variable = "Charlson_Comorbidity_Index",
                         subset = "admission",
                         test = "fisher"),
                    list(name = "RespSupp",
                         title = "Respiratory Support",
                         variable = "resp_supp_1.5",
                         subset = "all",
                         test = "fisher"),
                    list(name = "Revalidation",
                         title = "Days revalidation after discharge",
                         variable = "days_after_discharge_ICU_overall",
                         subset = "discharge",
                         test = "lmm"),
                    list(name = "AgeDischarge",
                         title = "Age",
                         variable = "age",
                         subset = "discharge",
                         test = "lmm"),
                    list(name = "AgeLowest",
                         title = "Age",
                         variable = "age",
                         subset = "lowest",
                         test = "lmm"),
                    list(name = "BMIAdmission",
                         title = "BMI",
                         variable = "BMI",
                         subset = "admission",
                         test = "lmm"),
                    list(name = "BMILowest",
                         title = "BMI",
                         variable = "BMI",
                         subset = "lowest",
                         test = "lmm"))

samples$`IL-8` <- samples$`IL-8_low_sens`

cytokines_APR <- c("IFN-gamma",
                   "TNF-alfa",
                   "IL-8",
                   "IL-6",
                   "IL-4",
                   "IL-2","IL-1Beta","IL-13","IL-12p70","IL-10",
                   "Eotaxin","Eotaxin_3","IP-10",
                   "MCP-1","MCP-4","MDC","MIP1a","MIP1b","TARC","GM-CSF",
                   "IL-12p40","IL-15","IL-16","IL-17","IL-1_alfa","IL-5",
                   "IL-7","TNFb","VEGF","IL1-RA","IL-18",
                   "CRP","Ddimers","Ferritin")
for(var in cytokines_APR){
  samples[[paste0("log10_", var)]] <- log10(as.numeric(samples[[var]]))
  to_evaluate[[length(to_evaluate)+1]] <-
    list(name = var,
         title = paste0("log10 ",var),
         variable = paste0("log10_", var),
         subset = "all",
         test = "lmm_randomEffect")
}


statistics <- data.frame(variable =character(),
                         profile1 = numeric(),
                         profile2 = numeric(),
                         pvalue = numeric(),
                         max = numeric())
for(evaluation in to_evaluate){
  subset <- subsets[[evaluation$subset]]$IDs
  statistics <- statistical_tests[[evaluation$test]](
    statistics, 
    name = evaluation$name,
    variable = samples[subset, evaluation$variable],
    clusters = samples[subset, "BloodProfile"],
    patientID = samples[subset, "PatientID"])
}
statistics <- finalize_statistics(statistics)
statistics_filtered <- filter_statistics(statistics)

plots <- list()
for(evaluation in to_evaluate){
  print(evaluation$name)
  p <- make_plot(title = evaluation$title,
                 variable = evaluation$variable,
                 subset = evaluation$subset)
  p <- add_stats(p, statistics_filtered, evaluation$name)
  
  ggsave(filename = paste0("Results/Figures/",
                           date, "_", clustername, "_", evaluation$name,".pdf"),
         width = 8, height = 6)
  
  plots[[evaluation$name]] <- p
}


plots[["Timepoint"]] <- plots[["Timepoint"]] + scale_y_discrete(expand = expansion(add = 1.5))
plots[["RespSupp"]] <- plots[["RespSupp"]] + scale_y_discrete(expand = expansion(add = 1.5))

pdf(paste0("Results/Figures/",
           date, "_", clustername, "_overviewStatistics.pdf"),
    width = 22, height = 12)
# Clinical
ggarrange(plotlist = plots[unique(statistics$variable)[!unique(statistics$variable) %in% cytokines_APR]])
# Significant cytokines
cytokines_to_plot <- sort(cytokines_APR[cytokines_APR %in% statistics_filtered$variable])

cytokines_to_plot[grep("IL", cytokines_to_plot)] <- cytokines_to_plot[grep("IL", cytokines_to_plot)][
  order(as.numeric(gsub("IL-([0-9]*).*", "\\1", cytokines_to_plot[grep("IL", cytokines_to_plot)])))]

gridExtra::grid.arrange(grobs = c(list(p_neutr_lymph = p_neutr_lymph),
                                     plots[cytokines_to_plot]),
                        layout_matrix = matrix(c(1,1, 2:5, 1, 1, 6:21), ncol = 6, byrow = TRUE))
dev.off()

xlsx::write.xlsx(x = statistics[,-5],
                 file = paste0("Results/Figures/",
                               date, "_", clustername, "_overviewStatistics.xlsx"))



# OTHER PLOTS ------------------------------------------------------------------

discharge <- samples[rownames(pctgs_highlevel), "rank"] == "7_samp_discharge"
p_revalidation_discharge_age <- 
  ggplot(samples[rownames(pctgs_highlevel),][discharge,],
         aes(x = as.numeric(age), 
             y = as.numeric(days_after_discharge_ICU_overall),
             col = clusters[discharge])) +
  ggbeeswarm::geom_quasirandom() +
  ggtitle("Days revalidation after discharge") + 
  xlab("Age (Discharge samples only)") + ylab("") +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  theme_minimal()


ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_RevalidationTime_Age.pdf"),
       width = 8, height = 6)


samples_unique <- 
  t(sapply(unique(as.character(samples$PatientID)),
           function(id) c(id,
                          samples[samples$PatientID == id, "max_resp_supp"][1],
                          samples[samples$PatientID == id, "BMI"][1],
                          samples[samples$PatientID == id, "age"][1],
                          samples[samples$PatientID == id, "gender"][1],
                          samples[samples$PatientID == id, "SOFA_first_24h"][1],
                          samples[samples$PatientID == id, "Charlson_Comorbidity_Index"][1]))) %>% 
  data.frame() %>% 
  magrittr::set_colnames(c("PatientID", "max_resp_supp", "BMI", 
                           "age", "gender", "SOFA_first_24h", "Charlson_Comorbidity_Index"))

for(i in c( "max_resp_supp", "BMI", "SOFA_first_24h", "Charlson_Comorbidity_Index")){
  samples_unique[[i]] <- as.numeric(samples_unique[[i]])
}
for(i in c( "max_resp_supp", "SOFA_first_24h", "Charlson_Comorbidity_Index")){
  samples_unique[[i]] <- factor(samples_unique[[i]])
}

samples_unique <- samples_unique[!is.na(samples_unique$max_resp_supp), ]


model <- lm(age ~ gender, 
            data = samples_unique)
emmeans <- emmeans(model, list(pairwise ~ gender), adjust = "none")
pdiff <- data.frame(emmeans$`pairwise differences of gender`)

p_age_gender <- ggplot(samples_unique,
                       aes(x = gender, 
                           y = as.numeric(age))) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  ggtitle(paste0("Gender / Age distribution")) +# (p = ", round(pdiff[,"p.value"], 2),")")) + 
  xlab("") + ylab("") +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  theme_minimal()


ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_age_gender.pdf"),
       width = 8, height = 6)


model <- lm(BMI ~ max_resp_supp, 
            data = samples_unique)
emmeans <- emmeans(model, list(pairwise ~ max_resp_supp), adjust = "none")
pdiff <- data.frame(emmeans$`pairwise differences of max_resp_supp`)

p_BMI_RespSupp <- ggplot(samples_unique,
                         aes(x = factor(max_resp_supp),
                             y = BMI)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  guides(col = FALSE) +
  xlab("Maximal Respiratory Support") + ylab("BMI") +
  ggtitle("BMI vs Maximal Respiratory Support")+
  theme_minimal()


ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_BMI_respSupp.pdf"),
       width = 8, height = 6)

model <- lm(BMI ~ SOFA_first_24h, 
            data = samples_unique)
emmeans <- emmeans(model, list(pairwise ~ SOFA_first_24h), adjust = "none")
pdiff <- data.frame(emmeans$`pairwise differences of SOFA_first_24h`)
pdiff <- cbind(pdiff,
               p_adjusted = p.adjust(pdiff[,"p.value"], method = "holm"))

p_BMI_SOFA <- ggplot(samples_unique,
                     aes(x = factor(SOFA_first_24h),
                         y = BMI)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  guides(col = FALSE) +
  xlab("SOFA first 24h") + ylab("BMI") +
  ggtitle("BMI vs SOFA")+
  theme_minimal()



ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_BMI_SOFA.pdf"),
       width = 8, height = 6)


comorbidity <- samples_unique$Charlson_Comorbidity_Index
max_resp_supp <- samples_unique$max_resp_supp
for(p1 in 1:4){
  for(p2 in (p1+1):5){
    fit <- fisher.test(comorbidity[max_resp_supp %in% c(p1,p2)], 
                       max_resp_supp[max_resp_supp %in% c(p1,p2)])
    print(paste0(p1, " - ", p2))
    print(fit$p.value)
    
  }
}
fisher.test(samples_unique$Charlson_Comorbidity_Index, 
            samples_unique$max_resp_supp)
p_Comorbidity_RespSupp <- ggplot(samples_unique,
                                 aes(x = factor(max_resp_supp),
                                     y = Charlson_Comorbidity_Index)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  guides(col = FALSE) +
  xlab("Maximal Respiratory Support") + ylab("Charlson Comorbidity Index") +
  ggtitle("Comorbidity vs Maximal Respiratory Support")+
  theme_minimal()


ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_Comorbity_respSupp.pdf"),
       width = 8, height = 6)


cor <- cor(as.numeric(samples[rownames(pctgs_highlevel),][discharge,"overall_stay_ICU"]),
           as.numeric(samples[rownames(pctgs_highlevel),][discharge,"days_after_discharge_ICU_overall"]),
           use = "complete")
p_los_discharge <- ggplot(samples[rownames(pctgs_highlevel),][discharge,],
                          aes(x = as.numeric(overall_stay_ICU),
                              y = as.numeric(days_after_discharge_ICU_overall),
                              col = clusters[discharge])) +
  geom_point() +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  xlab("ICU Days of stay (Discharge samples only)") + ylab("Days of revalidation") +
  ggtitle(paste0("ICU vs revalidation (r = ",round(cor, 3), ")"))+
  theme_minimal()


ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_LOS_revalidation.pdf"),
       width = 8, height = 6)





# Boxplot CD4/CD8 --------------------------------------------------------------

cd4 <- rowSums(pctgs_highlevel[,c("FlowSOM_CD4_naive",
                                  "FlowSOM_CD4_effector",
                                  "FlowSOM_CD4_Th1",
                                  "FlowSOM_CD4_Th2",
                                  "FlowSOM_CD4_Th17",
                                  "FlowSOM_CD4_Tfh",
                                  "FlowSOM_CD4_Treg")])

cd8 <- rowSums(pctgs_highlevel[,c("FlowSOM_CD8_naive",
                                  "FlowSOM_CD8_activated",
                                  "FlowSOM_CD8_memory")])


p_CD4CD8 <- ggplot(samples[rownames(pctgs_highlevel),],
                   aes(x = clusters,
                       y = cd4/cd8,
                       col = clusters)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  xlab("") + ylab("") +
  ggtitle("CD4/CD8 ratio")+
  theme_minimal()

ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_CD4CD8.pdf"),
       width = 8, height = 6)

# Days from symptom onset --------------------------------------------------------------


p_symptom <- ggplot(samples[rownames(pctgs_highlevel),],
                    aes(x = clusters,
                        y = as.numeric(Days_from_symp),
                        col = clusters)) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom() +
  guides(col = FALSE) +
  scale_color_manual(values = colors) +
  xlab("") + ylab("") +
  ggtitle("Days since onset of symptoms")+
  theme_minimal() +
  coord_flip()



ggsave(filename = paste0("Results/Figures/",
                         date, "_", clustername, "_symptoms.pdf"),
       width = 8, height = 6)




# Cytokines - Linear mixed model with random patient effect --------------------


pdf("Results/cytokine_normality.pdf", 
    width = 20, height = 20)
layout(matrix(1:72, ncol = 8, byrow = TRUE))
for(cytokine in cytokines){
  values <- as.numeric(samples[,cytokine])
  qqnorm(values, main = cytokine)
  qqline(values)
  
  qqnorm(log10(values), main = paste0("Log10 ",cytokine))
  qqline(log10(values))
}
dev.off()

cytokine_means <- list()
cytokine_statistics <- list()
for(cytokine in cytokines){
  print(cytokine)
  cytokine_data <- data.frame(Cytokine = as.numeric(samples[, cytokine]),
                              Log10Cytokine = log10(as.numeric(samples[, cytokine])),
                              BloodProfile = samples[, "BloodProfile"],
                              PatientID = samples[, "PatientID"])
  
  model <- lmer(Log10Cytokine ~ BloodProfile + (1|PatientID), 
                data = cytokine_data)
  emmeans <- emmeans(model, list(pairwise ~ BloodProfile), adjust = "none")
  cytokine_means[[cytokine]] <- data.frame(emmeans$`emmeans of BloodProfile`)[,"emmean"]
  pdiff <- data.frame(emmeans$`pairwise differences of BloodProfile`)
  df <- data.frame(cytokine = cytokine,
                   profile1 =  stringr::str_sub(pdiff[,"contrast"], 1, 1),
                   profile2 = stringr::str_sub(pdiff[,"contrast"], 5, 5),
                   pvalue = pdiff[,"p.value"])
  cytokine_statistics[[cytokine]] <- df
}
cytokine_means  <- do.call(rbind, cytokine_means)
cytokine_statistics <- do.call(rbind, cytokine_statistics)
cytokine_statistics$p_corrected <- p.adjust(cytokine_statistics$pvalue, 
                                            method = "holm")

cytokine_statistics$logfold_change <- 
  sapply(seq_len(nrow(cytokine_statistics)), function(i){
    cyt <- cytokine_statistics$cytokine[i]
    
    log10(cytokine_means[cyt, as.numeric(cytokine_statistics$profile1[i])] /
            cytokine_means[cyt, as.numeric(cytokine_statistics$profile2[i])] )
  })


# Volcano plot
dplyr::filter(cytokine_statistics, profile1 == 1) %>% 
  ggplot(aes(x = logfold_change,
             y = -log10(pvalue),
             label = cytokine,
             col = profile2)) +
  geom_text() +
  scale_color_manual(values = colors) +
  xlim(-0.5, 0.5)+
  theme_minimal()


cytokine_statistics$Stars <- ifelse(cytokine_statistics$p_corrected <= 0.05,
                                    ifelse(cytokine_statistics$p_corrected <= 0.01,
                                           ifelse(cytokine_statistics$p_corrected <= 0.001,
                                                  "***",
                                                  "**"),
                                           "*"),
                                    "")
filtered_p <- dplyr::filter(cytokine_statistics, p_corrected < 0.05)


plots <- list()
for(cytokine in cytokines){}
for(cytokine in unique(filtered_p$cytokine)){
  print(cytokine)
  values <- log10(as.numeric(samples[rownames(pctgs_highlevel), cytokine]))
  df <- data.frame(SampleID = rownames(pctgs_highlevel), 
                   Cytokine = values,
                   Cluster = clusters)
  plots[[cytokine]] <- ggplot(df) +
    geom_boxplot(aes(x = Cluster,
                       y = Cytokine,
                       col = Cluster),
                 outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom(aes(x = Cluster,
                                     y = Cytokine,
                                     col = Cluster)) +
    guides(col = FALSE) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(paste0("Log10 ",cytokine)) +
    xlab("") +
    ylab("")
  
  
  if(cytokine %in% filtered_p$cytokine){
    to_add <- filtered_p[filtered_p$cytokine == cytokine, ]
    to_add$y <- seq(from = max(df$Cytokine, na.rm = TRUE)*1.1,
                    by = max(df$Cytokine, na.rm = TRUE)*0.1,
                    length.out = nrow(to_add))
    plots[[cytokine]] <- plots[[cytokine]] +
      geom_segment(aes(x = profile1, xend = profile2,
                       y = y, yend = y), 
                   data = to_add) +
      geom_text(aes(x = (as.numeric(profile1)+as.numeric(profile2))/2,
                    y = 1.01*y,
                    label = Stars),
                data = to_add)
    
  }
  
  ggsave(filename = paste0("Results/Figures/",
                           date, "_", clustername, "_Cytokine_",cytokine,".pdf"),
         width = 8, height = 6)
}

ggarrange(plotlist = plots)



### Classical / non-classical monocytes ----------------------------------------

monocyte_df <- data.frame("Classical monocytes" = pctgs_highlevel[,"FlowSOM_classical_monocytes"],
                          "Non-classical monocytes" = pctgs_highlevel[,"FlowSOM_non-classical_monocytes"],
                          "SampleID" = rownames(pctgs_highlevel),
                          "BloodProfile" = samples[rownames(pctgs_highlevel), "BloodProfile"],
                          "rank" = samples[rownames(pctgs_highlevel), "rank"],
                          check.names = FALSE)
monocyte_df <- monocyte_df[order(clusters),]

pdf(file = paste0("Results/Figures/",
                      date, "_", clustername, "_RatioClassicalNonClassicalMonocytes.pdf"),
    width = 8, height = 6)
ggplot(monocyte_df) +
  geom_point(aes(x = `Classical monocytes`, 
                 y = `Non-classical monocytes`,
                 color = `BloodProfile`,
                 shape = rank == '11_healthy_control')) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, 
                                "FALSE" = 19)) + 
  guides(col = FALSE, shape = FALSE) +
  theme_minimal()


ggplot(monocyte_df) +
  geom_point(aes(x = 1:123, 
                 y = `Classical monocytes`/`Non-classical monocytes`, 
                 color = `BloodProfile`,
                 shape = rank == '11_healthy_control')) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, 
                                "FALSE" = 19)) + 
  guides(col = FALSE, shape = FALSE) +
  xlab("Samples ordered by blood profile 1-4") +
  theme_minimal()

p1 <- ggplot(monocyte_df) +
  geom_point(aes(x = 1:123, y = `Classical monocytes`, color = `BloodProfile`,
                 shape = rank == '11_healthy_control')) +
  guides(col = FALSE, shape = FALSE) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, 
                                "FALSE" = 19)) + 
  xlab("Samples ordered by blood profile 1-4") +
  theme_minimal()
p2 <- ggplot(monocyte_df) +
  geom_point(aes(x = 1:123, y = `Non-classical monocytes`, color = `BloodProfile`,
                 shape = rank == '11_healthy_control')) +
  guides(col = FALSE, shape = FALSE) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("TRUE" = 1, 
                                "FALSE" = 19)) + 
  xlab("Samples ordered by blood profile 1-4") +
  theme_minimal()

print(p1 | p2)

dev.off()

#### Overview change -------------------------------------------

samples$Date_of_sample_event <- strptime(samples$Date_of_sample_event, 
                                         format = "%d/%m/%Y__%H:%M:%S")

patients <- grep("HC", unique(samples$PatientID), value = TRUE, invert = TRUE)
blood_profiles <- list()
pattern_type <- list()
for(patient in patients){
  samples_patient <- samples[samples$PatientID == patient, ]
  samples_patient <- samples_patient[order(samples_patient$Date_of_sample_event), ]
  df <- cbind(samples_patient,
              BloodProfile = clusters[rownames(samples_patient)],
              Order = seq_len(nrow(samples_patient)))
  
  diff <-   diff(as.numeric(df$BloodProfile))
  
  pattern_type[[patient]] <- ifelse(all(diff == 0),
                                    "No change",
                                    ifelse(all(diff >= 0),
                                           "Improving",
                                           ifelse(all(diff <= 0),
                                                  "Decreasing",
                                                  "Variable")))
  df$pattern_type <-  pattern_type[[patient]]
  
  blood_profiles[[patient]] <- df
}

pattern_type <- unlist(pattern_type)
table(pattern_type)


blood_profiles_all <- do.call(rbind, blood_profiles)
rownames(blood_profiles_all) <- gsub(".*\\.", "", rownames(blood_profiles_all))

neutr_lymph_ratio <- pctgs_highlevel[rownames(blood_profiles_all), "FlowSOM_neutrophil"] / 
  pctgs_highlevel[rownames(blood_profiles_all), "FlowSOM_lymphocytes"]

p <- ggplot(blood_profiles_all,
            aes(x = as.numeric(Days_from_ad_ICU),
                y = BloodProfile)) +
  geom_line(aes(group = PatientID, col = pattern_type)) + 
  geom_point(aes(col = rank)) +
  theme_minimal() +
  facet_wrap(~ PatientID)



p <- ggplot(blood_profiles_all,
            aes(x = as.numeric(Days_from_ad_ICU),
                y = neutr_lymph_ratio)) +
  geom_line(aes(group = PatientID)) + 
  geom_point(aes(col = rank)) +
  theme_minimal() +
  facet_wrap(~ PatientID)


# Heatmap ----------------------------------------------------------------------

mfis <- GetMFIs(fsom_step2_lymph)
mfis <- mfis[,fsom_step2_lymph$FlowSOM$map$colsUsed]
colnames(mfis) <- GetMarkers(fsom_step2_lymph$FlowSOM, colnames(mfis))
rownames(mfis) <- paste("Cl",1:107,": ", cluster_label_2[1:107])

mfis_gr <- GetMFIs(fsom_step2_gran)
mfis_gr <- mfis_gr[,fsom_step2_gran$FlowSOM$map$colsUsed]
colnames(mfis_gr) <- GetMarkers(fsom_step2_gran$FlowSOM, colnames(mfis_gr))
rownames(mfis_gr) <- paste("grCl",1:100,": ", cluster_label_2[108:207])


pdf("Results/Figures/FlowSOM_heatmaps.pdf", width = 15, height = 30)
pheatmap::pheatmap(mfis)
pheatmap::pheatmap(mfis_gr)
dev.off()

# Old plots --------------------------------------------------------------------
blood_profile <- clus_ward

plots <- list()
for(clustering in list(pathsetter, clus_ward)){
  for(seed in 1:3){
    for(n_neighbors in 5:7){
      set.seed(seed)
      # umap <- uwot::umap(pctgs_highlevel, 
      #                    n_neighbors = n_neighbors,
      #                    min_dist = 0.5)
      # 
      # tsne <- Rtsne::Rtsne(pctgs_highlevel, 
      #                      perplexity = n_neighbors)
      
      pca <- prcomp(pctgs_highlevel)
      
      df <- data.frame(SampleID = rownames(pctgs_highlevel),
                       x = pca$x[,1], #tsne$Y[,1], #umap[,1],
                       y =pca$x[,2], #tsne$Y[,2],#umap[,2],
                       ratio = neutr_lymph_ratio,
                       blood_profile = clustering)
      plots[[length(plots) + 1]] <- 
        ggplot(df,
               aes(x = x,
                   y = y,
                   col = blood_profile)) +
        geom_point()  + 
        theme_minimal()#+
        #geom_text(aes(label = SampleID)) + 
        ggtitle(paste0("tSNE ",seed," ",n_neighbors))
    }
  }
}



pdf("Results/test_tsne.pdf", width = 12, height = 9)
ggarrange(plotlist = plots, ncol = 3, nrow = 3)
dev.off()


clus_umap <- cutree(hclust(dist(umap), method = "ward.D2"), k = 5)
# Order clusters by decreasing neutr/lymph ratio
clus_umap<- factor(as.numeric(factor(rank(-tapply(neutr_lymph_ratio, 
                                                  clus_umap,
                                                  mean))[as.character(clus_umap)])))

clusterings <-  list("Pathsetter" = pathsetter,
                     "FlowSOM" = clus_ward,
                     "UMAP" = clus_umap)

pdf("Results/200806_blood_profiles.pdf",
    width = 31, height = 15)
for(clustering in names(clusterings)){

  clusters <- clusterings[[clustering]]
  clusters[is.nan(clusters)] <- NA
  clusters <- factor(clusters)
  
  p_neutr_lymph <- ggplot(cbind(samples[rownames(pctgs_highlevel),],
                                pctgs_highlevel)[order(clusters),],
                          aes(x = 1:123, 
                              col = clusters[order(clusters)])) +
    geom_line(aes(y = FlowSOM_lymphocytes, group = 1), lwd = 1) +
    geom_line(aes(y = FlowSOM_neutrophil, group = 1), lwd = 1) +
    ggtitle("Neutr and Lymph") + 
    guides(col = FALSE) +
    xlab("Ordered by cluster") +
    ylab ("Neutrophils  / Lymphocytes") +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  
  
  discharge <- samples[rownames(pctgs_highlevel), "rank"] == "7_samp_discharge"
  p_revalidation_discharge <- 
    ggplot(samples[rownames(pctgs_highlevel),][discharge,],
           aes(x = clusters[discharge], 
               y = as.numeric(days_after_discharge_ICU_overall),
               col = clusters[discharge])) +
    geom_boxplot(outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom() +
    ggtitle("Overall revalidation discharge samples") + 
    guides(col = FALSE) +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  cluster_labels <- as.character(clusters)
  cluster_labels[is.na(cluster_labels)] <- "0"
  p_umap <- ggplot(samples[rownames(pctgs_highlevel),],
                   aes(x = umap[,1],
                       y = umap[,2],
                       col = clusters)) +
    geom_text(label = cluster_labels) + 
    ggtitle("UMAP") + 
    guides(col = FALSE) +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  pctgs_highlevel_long <- tidyr::gather(data.frame(SampleID = samples[rownames(pctgs_highlevel), "SampleID"],
                                                   Cluster = clusters,
                                                   pctgs_highlevel,
                                                   check.names = FALSE),
                                        "Celltype",
                                        "Pctgs",
                                        -SampleID, -Cluster) 
  pctgs_highlevel_long$Celltype <- factor(gsub("FlowSOM_", "", pctgs_highlevel_long$Celltype),
                                          levels = gsub("FlowSOM_", "", colnames(pctgs_highlevel)))
  p_celltypes <- ggplot(pctgs_highlevel_long, aes(x = Cluster,
                                         y = Pctgs,
                                         col = Cluster)) +
    geom_boxplot(outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap(~Celltype, scales = "free") +  
    guides(col = FALSE) +
    theme_minimal() +
    scale_color_manual(values = colors) +
    ggtitle(clustering)
  

  resp_supp <-  as.numeric(samples[rownames(pctgs_highlevel), "resp_supp_1.5"])
  resp_supp[is.na(resp_supp)] <- 0
  
  p_resp_supp <- ggplot(data.frame(SampleID = samples[rownames(pctgs_highlevel),"SampleID"], 
                                   RespSupp = resp_supp,
                                   Cluster = clusters),
                      aes(x = Cluster,
                          y = RespSupp,
                          col = Cluster,
                          group = Cluster)) +
    geom_boxplot(outlier.alpha = 0, col = "black") +
    ggbeeswarm::geom_quasirandom() +
    guides(col = FALSE) +
    ggtitle("Respiratory Support") +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  
  timepoint <- factor(samples[rownames(pctgs_highlevel),"rank"],
                      levels = unique(samples$rank)[
                        order(as.numeric(gsub("_.*", "", unique(samples$rank))))])
  p_timepoint <- ggplot(data.frame(SampleID = samples[rownames(pctgs_highlevel),"SampleID"], 
                                   Timepoint = timepoint,
                                   Cluster = clusters),
                        aes(x = Cluster,
                            y = Timepoint,
                            col = Cluster,
                            group = Cluster)) +
    geom_boxplot(outlier.alpha = 0, col = "black") +
    ggbeeswarm::geom_quasirandom() +
    guides(col = FALSE) +
    ggtitle("Timepoint of sample") +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  
  clinics_long_2 <- tidyr::gather(data.frame(SampleID = samples[rownames(pctgs_highlevel),"SampleID"], 
                                             CRP = as.numeric(samples[rownames(pctgs_highlevel),"CRP"]),
                                             Log_Ferritin = log(as.numeric(samples[rownames(pctgs_highlevel),"Ferritin"])),
                                             Log_Ddimers = log(as.numeric(samples[rownames(pctgs_highlevel),"Ddimers"])),
                                             Cluster = clusters),
                                  "Variable",
                                  "Value",
                                  -SampleID, -Cluster) 
  
  p_clinics2 <- ggplot(clinics_long_2, aes(x = Cluster,
                                           y = Value,
                                           col = Cluster)) +
    geom_boxplot(outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap(~Variable, scales = "free") +  
    guides(col = FALSE) +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  
  cytokines_long <- tidyr::gather(data.frame(SampleID = samples[rownames(pctgs_highlevel),"SampleID"], 
                                             logIL6 = log(as.numeric(samples[rownames(pctgs_highlevel),"IL-6"])),
                                             MCP1 = as.numeric(samples[rownames(pctgs_highlevel),"MCP-1"]),
                                             MCP4 = as.numeric(samples[rownames(pctgs_highlevel),"MCP-4"]),
                                             MDC = as.numeric(samples[rownames(pctgs_highlevel),"MDC"]),
                                             Cluster = clusters),
                                  "Variable",
                                  "Value",
                                  -SampleID, -Cluster) 
  
  p_cytokines <- ggplot(cytokines_long, aes(x = Cluster,
                                           y = Value,
                                           col = Cluster)) +
    geom_boxplot(outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap(~Variable, scales = "free", nrow = 1) +  
    guides(col = FALSE) +
    scale_color_manual(values = colors) +
    theme_minimal()
  
  annotation <- data.frame("BloodProfile" = factor(clusters))
  rownames(annotation) <- rownames(pctgs_highlevel)
  heatmap <- pheatmap::pheatmap(pctgs_highlevel[order(clusters), ],
                                gaps_row = cumsum(table(clusters)),
                                cluster_cols = FALSE,
                                cluster_rows = FALSE,
                                breaks = (1:100)^2 / (100)^2*100,
                                annotation_row = annotation,
                                annotation_colors = list(BloodProfile = colors),
                                silent = TRUE)
  
  print( ((p_neutr_lymph | p_umap) / ( p_revalidation_discharge | p_resp_supp | p_timepoint) / p_clinics2 / p_cytokines) | p_celltypes |  heatmap$gtable )
}

dev.off()


cytokines <- c("IFN-gamma",
               "TNF-alfa",
               "IL-8_low_sens",
               "IL-6",
               "IL-4",
               "IL-2","IL-1Beta","IL-13","IL-12p70","IL-10",
               "Eotaxin","Eotaxin_3","IL-8_high_sens","IP-10",
               "MCP-1","MCP-4","MDC","MIP1a","MIP1b","TARC","GM-CSF",
               "IL-12p40","IL-15","IL-16","IL-17","IL-1_alfa","IL-5",
               "IL-7","TNFb","VEGF","IL1-RA","IL-18")

log_cytokines <- c("IFN-gamma", "IL-6","IL-8_high_sens", "IL-18")
pdf("Results/cytokines.pdf",
    width = 31, height = 15)
for(clustering in names(clusterings)){
  
  clusters <- clusterings[[clustering]]
  clusters[is.nan(clusters)] <- NA
  clusters <- factor(clusters)
  plots <- list()
  for(cytokine in cytokines){
    values <- as.numeric(samples[rownames(pctgs_highlevel), cytokine])
    if(cytokine %in% log_cytokines) values <- log(values)
    df <- data.frame(SampleID = samples[rownames(pctgs_highlevel),"SampleID"], 
                     Cytokine = values,
                     Cluster = clusters)
    plots[[cytokine]] <- ggplot(df, aes(x = Cluster,
                                              y = Cytokine,
                                              col = Cluster)) +
      geom_boxplot(outlier.alpha = 0) +
      ggbeeswarm::geom_quasirandom() +
      guides(col = FALSE) +
      scale_color_manual(values = colors) +
      theme_minimal() +
      ggtitle(ifelse(cytokine %in% log_cytokines, paste0("log ",cytokine), cytokine))
  }
  
  print(ggarrange(plotlist = plots))
}

dev.off()


# Save to excel ----------------------------------------------------------------

names(clus_ward) <- rownames(pctgs_highlevel)
names(clus_umap) <- rownames(pctgs_highlevel)

samples_all <- readxl::read_xlsx("Metadata/Timeline_05082020_FDS.xlsx",
                                 sheet = "Final_table",
                                 na = "NA")
samples_all$FlowSOM_rank <- as.numeric(NA)
samples_all$FlowSOM_umap_rank <-  as.numeric(NA)
for(celltype in colnames(pctgs_highlevel)){
  samples_all[[celltype]] <-  as.numeric(NA)
}
for(sample_id in rownames(pctgs)){
  samples_all[which(samples_all$SampleID == sample_id), 
              c(paste0("Cl",1:107), paste0("grCl",1:100))] <-
    pctgs[sample_id, , drop = FALSE]
  samples_all[which(samples_all$SampleID == sample_id), 
              colnames(pctgs_highlevel)] <- pctgs_highlevel[sample_id, , drop = FALSE]
  samples_all[which(samples_all$SampleID == sample_id),
              "FlowSOM_rank"] <- as.numeric(clus_ward[sample_id])
  samples_all[which(samples_all$SampleID == sample_id),
              "FlowSOM_umap_rank"] <- as.numeric(clus_umap[sample_id])
}

xlsx::write.xlsx(samples_all,
                 "Metadata/Timeline_06082020_SVG.xlsx")
