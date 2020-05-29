library(ggpointdensity)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_2D_scatters <- function(fsom, 
                             channelpairs, 
                             ff,
                             max_bg_points = 3000, 
                             max_points = 1000,
                             size_bg_points = 0.5, 
                             size_points = 0.5,
                             metaclusters = NULL, 
                             clusters = NULL,
                             density = TRUE, 
                             centers = TRUE, 
                             color = NULL,
                             plot = TRUE, 
                             file_plot = "2D_scatter_plots.png", 
                             ...){
  
  i <- sample(nrow(fsom$FlowSOM$data), min(nrow(fsom$FlowSOM$data), max_bg_points))
  subsets <- c()
  if(length(metaclusters) > 0) subsets <- c(subsets, paste0("Metacluster ", metaclusters))
  if(length(clusters) > 0)  subsets <- c(subsets, paste0("Cluster ", clusters)) # loop over all subsets at once
  plots_list <- list()
  for (subset in subsets){
    for (channelpair in channelpairs){
      # background dataframe
      df_bg <- data.frame(fsom$FlowSOM$data[i, c(channelpair[1], channelpair[2])]) 
      colnames(df_bg) <- c("m1", "m2")
      
      # (meta)cluster dataframe
      n <- as.numeric(sub(".* ", "", subset))
      if (startsWith(subset, "Meta")){ # metaclusters
        df_ss <- fsom$FlowSOM$data[which(fsom$metaclustering[fsom$FlowSOM$map$mapping[, 1]]==n), c(channelpair[1], channelpair[2])]
        df_c <- data.frame(fsom$FlowSOM$map$medianValues[which(fsom$metaclustering[1:nrow(fsom$FlowSOM$map$medianValues)]==n), 
                                                         c(channelpair[1], channelpair[2]),
                                                         drop = FALSE])
        col <- gg_color_hue(length(levels(fsom$metaclustering)))[n]
      } else { # clusters
        df_ss <- fsom$FlowSOM$data[fsom$FlowSOM$map$mapping[, 1]==n, c(channelpair[1], channelpair[2])]
        df_c <- data.frame(fsom$FlowSOM$map$medianValues[n, 
                                                         c(channelpair[1], channelpair[2]),
                                                         drop = F])
        col <- gg_color_hue(length(levels(fsom$metaclustering)))[fsom$metaclustering[n]]
      }
      df_ss <- data.frame(df_ss[sample(nrow(df_ss), min(nrow(df_ss), max_points)),])
      colnames(df_ss) <- c("m1", "m2")
      colnames(df_c) <- c("m1", "m2")
      
      p <- ggplot(data = df_ss, aes(x = m1, y = m2)) +
        geom_point(data = df_bg, colour = "grey", size = size_bg_points) + # background dot plot
        theme_classic() +
        ggtitle(subset) +
        xlab(get_markers(ff, channelpair[1])) +
        ylab(get_markers(ff, channelpair[2])) + 
        theme(legend.position = "none")
      
      if (density) {
        p <- p  + geom_pointdensity(adjust = 20, size = size_points) + # subset density plot
          scale_color_gradientn(colors = colorRamps::matlab.like2(10))
      } else {
        if (is.null(color)){
          p <- p + geom_point(colour = col, size = size_points) #subset plot, metacluster colors
        } else {
          p <- p + geom_point(colour = color[match(subset, subsets)], size = size_points) #subset plot, metacluster colors
        }
        
      }
      
      if (centers) { # cluster centers
        p <- p + geom_point(data = df_c, shape=21, fill="white", color="black", size=3)
      }
      
      plots_list[[length(plots_list)+1]] <- p
    }
  }
  
  if (plot) {
    png(file_plot, width = 400 * length(channelpairs), height = 400 * length(subsets))
    print(ggarrange(plotlist = plots_list, 
                    ncol = length(channelpairs), nrow = length(subsets), 
                    common.legend = F))
    dev.off()
  } else {
    return(plots_list)
  }
}