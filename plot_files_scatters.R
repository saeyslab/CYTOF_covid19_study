library(ggplot2)
library(ggpubr)

plot_files_scatters <- function(input, 
                                channels = NULL, 
                                y_margin = NULL, 
                                names = NULL,
                                labels = NULL, 
                                color = NULL, 
                                legend = FALSE,
                                max_points = 50000,
                                ncol = NULL, 
                                nrow = NULL,
                                plot = TRUE, 
                                file_plot = "files_scatters.png"){
  
  if (is(input, "flowSet")) {
    data <- flowCore::fsApply(input, function(ff) {
      ff@exprs
    })
    cell_counts <- flowCore::fsApply(input, function(ff) {
      nrow(ff)
    })
    file_values <- unlist(sapply(seq_len(length(cell_counts)), 
                                 function(i) {
                                   rep(i, cell_counts[i])
                                 }))
    ff <- input[[1]]
  } else if (is(input, "flowFrame")) {
    ff <- input
    data <- ff@exprs
    file_values <- data[, "File"]
  } else {
    ff <- AggregateFlowFrames(input,
                              cTotal = max_points)
    data <- ff@exprs
    file_values <- data[, "File"]
  }
  n_files  <- length(unique(file_values))
  
  subset <- sample(seq_len(nrow(data)), min(max_points, nrow(data)))
  if (is.null(channels)) {
    data <- data[subset,] } else {
      data <- data[subset, channels]
    }
  file_values <- file_values[subset]
  channels <- colnames(data)
  
  if (is.null(names)) { # if no names are provided, the files will be numbered
    names <- as.character(seq_len(n_files))
  }
  
  if (is.null(labels)) { # if there are no groups, all files will be labeled "1"
    labels <- rep("1", n_files)
  }
  
  plots_list <- list()
  for (channel in channels) {
    df <- data.frame("intensity" = data[,channel],
                     "names" = factor(names[file_values], levels = unique(names)),
                     "label" = factor(labels[file_values], levels = unique(labels)))
    p <- ggplot(df, aes(names, intensity)) +
      geom_jitter(position = position_jitter(width=0.1), alpha = 0.5, aes(colour = label), shape = ".") +
      ylab(FlowSOM::get_markers(ff, channel)) +
      theme_classic() +
      theme(axis.text.x=element_text(angle =- 90, hjust = 0, v = 0.5),
            axis.title.x=element_blank()) + 
      guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1)))
    
    if (!is.null(color)) { # if manual colors are provided
      p <- p + scale_color_manual(values = color)
    }
    
    if (!is.null(y_margin)) { # if y margins are provided
      p <- p + ylim(y_margin)
    }
    
    if (!legend) { # if you don't want a legend on the plot
      p <- p + theme(legend.position = "none")
    }
    
    plots_list[[length(plots_list)+1]] <- p
  }
  
  if (plot) {
    if (is.null(nrow) & is.null(ncol)) {
      nrow <- floor(sqrt(length(channels)))
      ncol <- ceiling(length(channels)/nrow)
    } else if (!is.null(nrow) & !is.null(ncol)) {
      if(nrow*ncol < length(channels)) (stop("too few rows/cols to make plot"))
    } else if (is.null(nrow)) {
      nrow <- ceiling(length(channels)/ncol)
    } else {
      ncol <- ceiling(length(channels)/nrow)
    }
    
    png(file_plot, width = ncol*(15+15*n_files), height = 300*nrow)
    p <- annotate_figure(ggarrange(plotlist = plots_list,
                                   common.legend = legend, 
                                   ncol = ncol, nrow = nrow),
                         bottom = text_grob("Files"))
    print(p)  
    dev.off()
  } else {
    return(plots_list)
  }
}