
add_count <- function(df, column){
  count <- table(df[, column])
  df[, paste0(column, "_count")] <- paste0(df[, column], 
                                           " (", count[df[, column]], ")")
  if(is.factor(df[, column])) {
    df[, paste0(column, "_count")] <- 
      factor(df[, paste0(column, "_count")],
             levels = paste0(levels(df[, column]),
                             " (", count[levels(df[, column])], ")"))
  }
  return(df)
}

get_markers_per_file <- function(samples,
                                 result_rds = "RDS/markers_per_file.RDS",
                                 recompute = FALSE){
  if (!recompute & file.exists(result_rds)) {
    markers_per_file <- readRDS(result_rds)
    message("Reloaded previous result for markers_per_file. ",
            "Set recompute = TRUE to recompute." )
  } else {
    markers_per_file <- list()
    for (sample in rownames(samples)) {
      ff <- flowCore::read.FCS(samples[sample, "File"],
                     emptyValue = FALSE)
      markers_per_file[[sample]] <- FlowSOM::GetMarkers(ff, 
                                                        flowCore::colnames(ff))
    }
    saveRDS(markers_per_file,
            result_rds)
  }
  
  return(markers_per_file)
}

load_marker_excel <- function(excel_file = "Metadata/markers_immune_cytof_panel_v4 (1).xlsx",
                              sheet_name = "Extended panel"){
  
  marker_info <- xlsx::read.xlsx(excel_file,
                                 sheetName = sheet_name,
                                 check.names = FALSE,
                                 stringsAsFactors = FALSE)
  
  marker_info$RefName <- paste0(marker_info$Isotope, "_", marker_info$Marker)
  
  alternative_names <- c("141Pr_CCR6" = "141Pr_CD196_CCR6",
                         "143Nd_CD123" = "143Nd_CD123_IL-3R",
                         "146Nd_CD8" = "146Nd_CD8a",
                         "152Sm_CCR4" = "152Sm_CD194_CCR4",
                         "153Eu_CD25" = "153Eu_CD25_IL-2Ra",
                         "156Gd_CXCR3" = "156Gd_CD183_CXCR3",
                         "158Gd_CXCR5" = "158Gd_CD185_CXCR5",
                         "163Dy_CD56" = "163Dy_CD56_NCAM",
                         "164Dy_TCRyd" = "164Dy_TCRgd",
                         "166Er_CD294/CRTH2" = "166Er_CD294",
                         "167Er_CCR7" = "167Er_CD197_CCR7",
                         "169Tm_CD159a/NKG2A" = "169Tm_NKG2A",
                         "175Lu_CD279 [PD-1]" = "175Lu_PD-1",
                         "176Yb_CD127" = "176Yb_CD127_IL-7Ra")
  
  for(m in names(alternative_names)){
    marker_info$RefName[marker_info$RefName == m] <- alternative_names[m]
  }
  
  marker_info$Channel <- gsub("([0-9]*)([A-z]*)", "\\2\\1Di", marker_info$Isotope)
  marker_info$Channel[marker_info$Channel == "TBADi"] <- NA
  
  rownames(marker_info) <- marker_info$RefName
  return(marker_info)
}

library(dendextend) 

plot_cluster_dendrogram <- function(hierarchical_clustering,
                                    values,
                                    values_title,
                                    leaf_col,
                                    main){
  
  values_col <- rev(RColorBrewer::brewer.pal(8, "RdYlBu"))[
    cut(values, 
        breaks = 7,
        labels = FALSE)]
  
  
  get_h_properties <- function(x){
    if(!is.null(attr(x, "leaf"))) {
      return( list(height = attr(x, "height")))
    } else {
      return( list(height = attr(x, "height"),
                   child1 = get_h_properties(x[[1]]),
                   child2 = get_h_properties(x[[2]])))
    }
  }
  dend <- as.dendrogram(hierarchical_clustering)
  
  heights <- get_h_properties(dend)
  heights <- unlist(heights)
  heights_o <- order(heights)
  nodes_col <- rep("black", length(values)*2 + 1)
  nodes_col[heights_o[(length(values)+2):(length(values)*2+1)]] <- values_col[seq_len(length(values))]
  
  layout(matrix(c(1,rep(2, 3),
                  3, rep(2, 3)), nrow = 2, byrow = TRUE))
  plot(values,
       col = values_col,
       pch = 19,
       ylab = values_title,
       xlab = "",
       xaxt = "n")
  
  plot(dend %>% 
         set("nodes_pch", 19) %>%   
         set("nodes_col", nodes_col) %>% 
         set("labels_cex", 0.5) %>%
         set("leaves_col", leaf_col[hierarchical_clustering$order]), 
       horiz = TRUE,
       main = main)
  
}

compute_quantiles <- function(samples,
                              markers,
                              aggregate,
                              recompute = FALSE){
  quantile_file <- "RDS/quantiles.RDS"
  if(!recompute & file.exists(quantile_file)){
    quantiles <- readRDS(quantile_file)
    message("Reusing quantiles. Set recompute = TRUE to recompute.")
  } else {
    
    quantiles <- matrix(NA, 
                        nrow = nrow(samples),
                        ncol = 0,
                        dimnames = list(samples$Sample_ID,
                                        NULL))
    for(channel in names(markers)){
      x <- aggregate@exprs[, channel]
      q <- tapply(x, aggregate@exprs[,"File"],
                  quantile, c(0.25, 0.5, 0.75))
      q <- do.call(rbind, q)
      colnames(q) <- paste(c("q25", "q50", "q75"), markers[channel])
      quantiles <- cbind(quantiles, q)
    }
    
    quantiles <- quantiles[, apply(quantiles, 2, var) != 0]
    saveRDS(quantiles, quantile_file)
  }
  return(quantiles)
}

plot_quantiles <- function(quantiles){
  
  quantile_annotation <- data.frame(q = gsub(" .*", "", colnames(quantiles)),
                                    marker = gsub(".*_", "", colnames(quantiles)))
  rownames(quantile_annotation) <- colnames(quantiles)
  quantile_annotation <- ComplexHeatmap::HeatmapAnnotation(df = quantile_annotation)
  
  p <- ComplexHeatmap::Heatmap(abs(scale(quantiles)),
                               col = c("white", "#1d91c0"),
                               cluster_columns = FALSE,
                               column_split = gsub("_", "\n", gsub(".* ", "", colnames(quantiles))),
                               column_title_gp = grid::gpar(fontsize = 6),
                               column_names_gp = grid::gpar(fontsize = 6),
                               row_names_gp = grid::gpar(fontsize = 7),
                               column_labels = gsub(" .*", "", colnames(quantiles)))
  
  pca <- prcomp(scale(quantiles))
  p2 <- ggplot(data.frame(pca1 = pca$x[,1],
                          pca2 = pca$x[,2],
                          Sample_ID = rownames(quantiles))) +
    geom_point(aes(x = pca1, y = pca2), col = "#1d91c0") +
    geom_text(aes(x = pca1, y = pca2, label = Sample_ID), size = 3) +
    theme_minimal()
  
  return(list(p, p2))
}



# Stats ------------------------------------------------------------------------
# 
# add_fisher <- function(statistics,
#                        name,
#                        variable,
#                        clusters,
#                        patientID){
#   for(p1 in 1:3){
#     for(p2 in (p1+1):4){
#       fit <- fisher.test(variable[clusters %in% paste0("BP", c(p1,p2))], 
#                          clusters[clusters %in% paste0("BP", c(p1,p2))])
#       
#       statistics <- rbind(statistics,
#                           data.frame(variable = name,
#                                      profile1 = paste0("BP", p1),
#                                      profile2 = paste0("BP", p2),
#                                      pvalue = fit$p.value,
#                                      max = max(as.numeric(variable), na.rm = TRUE),
#                                      range = diff(range(as.numeric(variable), na.rm = TRUE))))
#       
#     }
#   }
#   return(statistics)
# }

# library(brms)
# library(tidybayes)
# add_brm <- function(statistics, 
#                                  name,
#                                  variable,
#                                  clusters,
#                                  patientID){
#   
#   data <- data.frame(Variable = variable,
#                      BloodProfile = clusters,
#                      PatientID = patientID)
#   
#   model <- brm(Variable ~ BloodProfile + (1|PatientID), family = "categorical",
#                  data = data, prior = c(set_prior("normal (0, 8)")), 
#                control = list(adapt_delta = 0.99))
#   
#   
#   
#   model <- brm(Variable ~ BloodProfile + (1|PatientID), 
#                family = "categorical",
#                data = data, 
#                prior = set_prior("normal (0, 1)"))
#   # Problem with post hoc test
#   
#   emmeans <- emmeans(model, pairwise ~ BloodProfile, 
#                      adjust = "none", dpar = "mu4sampadICU") 
#   cont <- contrast(emmeans)
#   
#   cont_posterior <- gather_emmeans_draws(cont)
#   
#   
#   emmeans <- emmeans(model, pairwise ~ BloodProfile, adjust = "none", dpar = "muHFNC") 
#   # 
#   # pdiff <- data.frame(emmeans$`pairwise differences of BloodProfile`)
#   # statistics <- rbind(statistics,
#   #                     data.frame(variable = name,
#   #                                profile1 =  stringr::str_sub(pdiff[,"contrast"], 1, 3),
#   #                                profile2 = stringr::str_sub(pdiff[,"contrast"], 7, 9),
#   #                                pvalue = pdiff[,"p.value"],
#   #                                max = max(variable, na.rm = TRUE),
#   #                                range = diff(range(variable, na.rm = TRUE))))
#   # return(statistics)
# }

add_lmm <- function(statistics,
                    name,
                    variable,
                    clusters,
                    patientID){
  
  data <- data.frame(Variable = variable,
                     BloodProfile = clusters)
  model <- lm(Variable ~ BloodProfile, 
              data = data)
  emmeans <- emmeans(model, list(pairwise ~ BloodProfile), adjust = "none")
  pdiff <- data.frame(emmeans$`pairwise differences of BloodProfile`)
  statistics <- rbind(statistics,
                      data.frame(variable = name,
                                 profile1 =  stringr::str_sub(pdiff[,"X1"], 1, 2),
                                 profile2 = stringr::str_sub(pdiff[,"X1"], 6, 7),
                                 pvalue = pdiff[,"p.value"],
                                 max = max(variable, na.rm = TRUE),
                                 range = diff(range(variable, na.rm = TRUE)),
                                 method = "lmm"))
  return(statistics)
}

add_lmm_randomEffect <- function(statistics,
                                 name,
                                 variable,
                                 clusters,
                                 patientID){
  
  data <- data.frame(Variable = variable,
                     BloodProfile = clusters,
                     PatientID = patientID)
  
  model <- lmer(Variable ~ BloodProfile + (1|PatientID), 
                data = data)
  emmeans <- emmeans(model, list(pairwise ~ BloodProfile), adjust = "none")
  
  pdiff <- data.frame(emmeans$`pairwise differences of BloodProfile`)
  statistics <- rbind(statistics,
                      data.frame(variable = name,
                                 profile1 =  stringr::str_sub(pdiff[,"X1"], 1, 2),
                                 profile2 = stringr::str_sub(pdiff[,"X1"], 6, 7),
                                 pvalue = pdiff[,"p.value"],
                                 max = max(variable, na.rm = TRUE),
                                 range = diff(range(variable, na.rm = TRUE)),
                                 method = "lmm_randomEffect"))
  return(statistics)
}

add_lmm_randomEffect_means <- function(means,
                                       name,
                                       variable,
                                       clusters,
                                       patientID){
  
  data <- data.frame(Variable = variable,
                     BloodProfile = clusters,
                     PatientID = patientID)
  
  model <- lmer(Variable ~ BloodProfile + (1|PatientID), 
                data = data)
  emmeans <- emmeans(model, list(pairwise ~ BloodProfile), adjust = "none")
  emmeans  <- data.frame(emmeans$`emmeans of BloodProfile`)
  #pdiff <- data.frame(emmeans$`pairwise differences of BloodProfile`)
  means <- rbind(means,
                 data.frame(variable = name,
                            profile1 = emmeans[1, "emmean"],
                            profile2 = emmeans[2, "emmean"],
                            profile3 = emmeans[3, "emmean"],
                            profile4 = emmeans[4, "emmean"]))
  return(means)
}
statistical_tests <- list("lmm" = add_lmm,
                          "lmm_randomEffect" = add_lmm_randomEffect)


finalize_statistics <- function(statistics){
  statistics$p_corrected <- p.adjust(statistics$pvalue, 
                                     method = "BH")
  
  statistics$stars <- ifelse(statistics$p_corrected <= 0.05,
                             ifelse(statistics$p_corrected <= 0.01,
                                    ifelse(statistics$p_corrected <= 0.001,
                                           "***",
                                           "**"),
                                    "*"),
                             "")
  return(statistics)
}

filter_statistics <- function(statistics){
  statistics_filtered <- dplyr::filter(statistics, p_corrected < 0.05)
  
  statistics_filtered$y <- make.unique(statistics_filtered$variable) %>% 
    gsub(".*\\.([0-9]*)", "\\1", .) %>% 
    as.numeric() %>% 
    (function(x){x[is.na(x)] <- 0; x})() %>% 
    (function(x){statistics_filtered$max + (1+x) * statistics_filtered$range*0.05})()
  
  return(statistics_filtered)
}

make_plot <- function(title,
                      variable,
                      subset){
  subset_ids <- subsets[[subset]]$IDs
  ggplot(samples[subset_ids, ]) +
    geom_boxplot(aes_string(x = "BloodProfile",
                            y = paste0("`", variable, "`"),
                            col = "BloodProfile",
                            group = "BloodProfile"),
                 outlier.alpha = 0) +
    ggbeeswarm::geom_quasirandom(aes_string(x = "BloodProfile",
                                            y = paste0("`", variable, "`"),
                                            col = "ColorType",
                                            shape = "rank == 'Healthy control'"))+
                                 #size = 3) +
    guides(col = FALSE, shape = FALSE) +
    ggtitle(title) +
    xlab(subsets[[subset]]$Description) + ylab("") +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c("TRUE" = 1, 
                                  "FALSE" = 19)) + 
    theme_minimal() +
    theme(plot.margin = margin(t = 5)) +
    coord_cartesian(clip = "off")
}

add_stats <- function(p, statistics, stat_name){
  if(stat_name %in% statistics$variable){
    p <- p +
      geom_segment(aes(x = profile1, xend = profile2,
                       y = y,  yend = y), 
                   data = statistics[statistics$variable == stat_name, , drop = FALSE]) +
      geom_text(aes(x = ((5-as.numeric(gsub("R", "", profile1)))+(5-as.numeric(gsub("R", "", profile2))))/2,
                    y = 1.01*y,
                    label = stars),
                data = statistics[statistics$variable == stat_name, , drop = FALSE]) 
  }
  return(p)
}
