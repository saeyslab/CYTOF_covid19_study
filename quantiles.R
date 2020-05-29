library(ComplexHeatmap)

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
                               column_title_gp = gpar(fontsize = 6),
                               column_names_gp = gpar(fontsize = 6),
                               row_names_gp = gpar(fontsize = 7),
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