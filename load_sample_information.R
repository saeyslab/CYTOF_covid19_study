library(tidyverse)
library(flowCore)

#' Load Sample Information
#' 
#' List all fcs files in a directory and parse their patient id and timepoint
#' @param meta_file Excel file with meta information
#' @param fcs_dir Path to the directory to use. Default "All_cleaned_fcs_files"
#' @return Dataframe with four columns: File, Patient_ID, Timepoint and 
#'         Sample_ID (format "Patient_ID Timepoint", 
#'         to which a dot and a number is added to make unique if necessary). 
#'         Rownames correspond to Sample_ID.
#' @examples
#' samples <- load_sample_information(fcs_dir = "PeacoQC_results/fcs_files",
#'                                    plot = "SampleOverview_PeacoQC.pdf")

load_sample_information <- function(meta_file,
                                    selection_file,
                                    fcs_dir,
                                    recompute = FALSE) {
  sample_file <- "RDS/samples.RDS"
  if (!recompute & file.exists(sample_file)) {
    samples <- readRDS(sample_file)
    message("Reloaded previous result. Use recompute = TRUE to recompute.")
    return(samples)
  } else {
    samples <- xlsx::read.xlsx(meta_file,
                               sheetName = "Sheet1",
                               check.names = FALSE,
                               stringsAsFactors = FALSE)
    samples <- distinct(samples)
    
    samples_of_interest <- xlsx::read.xlsx(selection_file,
                                           sheetName = "Groups",
                                           check.names = FALSE,
                                           stringsAsFactors = FALSE)
    samples_of_interest <- samples_of_interest[grep(1, samples_of_interest$`measured_in_cytof?`),]
    
    samples <- dplyr::left_join(x = samples,
                                y = samples_of_interest[,1:11],
                                by = c("File_Names" = "cleaned_file_name"))
    
    samples$File_peacoQC <- NA
    samples$nCells <- NA
    samples$Date_measured <- NA
    
    for (i in seq_len(nrow(samples))){
      file <- samples$File_Names[i]
      file_peacoQC <- file.path(fcs_dir,
                                paste0(sub(".fcs", "", basename(file)), 
                                       "_QC.fcs"))
      if(file.exists(file_peacoQC)){
        samples[i, "File_peacoQC"] <- file_peacoQC
        ff <- read.FCS(file_peacoQC)
        samples[i, "nCells"] <- nrow(ff)
        samples[i, "Date_measured"] <- paste(ff@description$`$DATE`, ff@description$`$BTIM`)
        #format = "%d-%b-%Y %H:%M:%S"
      }
    }
    
    samples$Sample_ID <- paste0(samples$patient_ID, " ", samples$subtype)
    samples$Sample_ID <- make.unique(samples$Sample_ID)
    
    samples <- samples[order(samples$patient_ID, samples$subtype), ]
    
    rownames(samples) <- samples$Sample_ID
    saveRDS(samples,
            sample_file)
  }
  
  return(samples)
}

plotSampleOverview <- function(samples,
                               plot = NULL){
  
  subtype_count <- table(samples$subtype[!(is.na(samples$File_peacoQC) | samples$nCells < 50000)])
  samples$subtype_count <- paste0(samples$subtype, 
                                  " (", subtype_count[samples$subtype],")")
  
  patient_count <- table(samples$patient_ID[!(is.na(samples$File_peacoQC) | samples$nCells < 50000)])
  samples$patient_ID_count <- paste0(samples$patient_ID, 
                                     " (", patient_count[samples$patient_ID],")")
  samples$patient_ID_count <- factor(as.character(samples$patient_ID_count), 
                                     levels = unique(samples$patient_ID_count[order(samples$subtype,
                                                                                    samples$patient_ID)]))

  
  samples$used_to_train <- 
    as.character(! is.na(samples$File_peacoQC) & 
                   samples$nCells >= 50000 &
                 ! samples$Sample_ID %in%
                   c("COVID_ICU_013 4_ICU_A",
                     "COVID_ICU_014 4_ICU_A",
                     "COVID_ICU_015 4_ICU_A",
                     "COVID_ICU_016 4_ICU_A",
                     "COVID_ICU_003 8_ICU_D",
                     "COVID_W_020 1_W_O",
                     "COVID_W_023 1_W_O",
                     "COVID_HC_002 0_CTR",
                     grep("\\.1", samples$Sample_ID, value = TRUE),
                     grep("HEALTHY", samples$Sample_ID, value = TRUE)))
  
  p <- ggplot(samples) +
    geom_point(aes(x = patient_ID_count, 
                   y = subtype_count, 
                   col = #ifelse(is.na(File_peacoQC) | nCells < 50000 | (!is.na(QC_cytof) & QC_cytof == "0"),
                          #      "grey",
                                ifelse(is.na(Clinical_condition_at_sampling),
                                       "black",
                                       ifelse(Clinical_condition_at_sampling == "good",
                                              "#7fcdbb",
                                              ifelse(Clinical_condition_at_sampling == "bad",
                                                     "#d95f0e", 
                                                     "#dd1c77"))),  # "healthy"
                   size = used_to_train)) +
    scale_color_identity() +
    scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -90, v = 0.5, hjust = 0, size = 7),
          legend.position = "none")
  
  if (!is.null(plot)){
    ggsave(filename = plot, 
           plot = p,
           width = 11.7,
           height = 4)
  }
  
  return(p)
}