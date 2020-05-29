source("R_scripts/plot_files_scatters.R")

make_aggregate <- function(samples,
                           channels,
                           cTotal = 3000000,
                           agg_file = "Aggregate/aggregate.fcs",
                           seed = NULL,
                           recompute = FALSE){
  
  if (!dir.exists(dirname(agg_file))) 
    dir.create(dirname(agg_file), recursive = TRUE)
  
  if (!recompute & file.exists(agg_file)) {
    
    aggregate <- read.FCS(agg_file)
    message("Reloaded previously generated aggregate.\n",
            "Set recompute = TRUE or choose a different agg_file to recompute.")
    
  } else {
    
    if (!is.null(seed)) set.seed(seed)
    aggregate <- FlowSOM::AggregateFlowFrames(fileNames = samples$File_peacoQC,
                                              cTotal = cTotal,
                                              channels = channels,
                                              writeOutput = TRUE,
                                              outputFile = agg_file,
                                              writeMeta = TRUE,
                                              verbose = TRUE)
  }
  
  return(aggregate)
}