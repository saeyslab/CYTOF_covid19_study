#install.packages("dendextend")
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