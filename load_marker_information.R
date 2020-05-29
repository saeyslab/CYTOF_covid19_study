load_marker_excel <- function(excel_file = "Meta_information/markers_immune_cytof_panel_v4 (1).xlsx",
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
      if(!is.na(samples[sample, "File_peacoQC"])){
        ff <- read.FCS(samples[sample, "File_peacoQC"],
                       emptyValue = FALSE)
        markers_per_file[[sample]] <- FlowSOM::get_markers(ff, colnames(ff))
      }
    }
    saveRDS(markers_per_file,
            result_rds)
  }
  
  return(markers_per_file)
}

plotMarkerOverview <- function(samples,
                               markers_measured,
                               plot = NULL){
  
  sample_order <- order(as.POSIXct(samples$Date_measured,
                                   format = "%d-%b-%Y %H:%M:%S"))
  markers_measured$Sample_ordered <- factor(markers_measured$Sample,
                                            levels = samples$Sample_ID[sample_order])
  sample_colors <- ifelse(samples$subtype[sample_order] %in% c("1_W_O", "4_ICU_A"),
                          "#1d91c0", "black")
  
  p <- ggplot(markers_measured) +
    geom_point(aes(x = Sample_ordered, 
                   y = potentialMarker, 
                   col = ifelse(!channel_of_interest,
                                "grey",
                                ifelse(channel_type == "DIP panel", 
                                       "#1d91c0", "black"))), 
               size = 2) +
    scale_color_identity() +
    scale_y_discrete(limits = rev(levels(markers_measured$potentialMarker))) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -90, v = 0.5, hjust = 0,
                                     color = sample_colors,
                                     size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
  
  if(!is.null(plot))
    ggsave(filename = plot,
           plot = p,
           width = 11.7, height = 8.3)
  
  return(p)
}

load_markers_measured <- function(samples,
                                  excel_file,
                                  sheet_name,
                                  recompute){
  
  marker_info <- load_marker_excel(excel_file, sheet_name)
  
  markers_per_file <- get_markers_per_file(samples,
                                           result_rds = "RDS/markers_per_file.RDS",
                                           recompute)
  
  markers_measured <- lapply(names(markers_per_file), 
                             function(x) data.frame(Sample = x,
                                                    Channel = names(markers_per_file[[x]]),
                                                    Marker = markers_per_file[[x]],
                                                    stringsAsFactors = FALSE))
  markers_measured <- do.call(rbind, markers_measured)
  rownames(markers_measured) <- NULL
  markers_measured$Marker[markers_measured$Marker == "PD-1"] <- "175Lu_PD-1"
  
  markers_measured$channel_of_interest <- 
    grepl("_", markers_measured$Marker) & ! grepl("Bead|Event|Eq4|DNA|Live", markers_measured$Marker)
  markers_measured$channel_type <- factor(sapply(markers_measured$Channel,
                                                 function(x){
                                                   ifelse(x %in% marker_info$Channel, 
                                                          marker_info$`CyTOF panel`[grep(x, marker_info$Channel)], 
                                                          "Other")}),
                                          levels = c("DIP panel", "added", "Other"))
  markers_measured$potentialMarker <- sapply(markers_measured$Channel, function(x) paste(unique(markers_measured$Marker[markers_measured$Channel == x]), collapse = " / "))
  markers_measured$potentialMarker <- factor(markers_measured$potentialMarker,
                                             levels = unique(markers_measured$potentialMarker[order(markers_measured$channel_type)]))
  
  return(markers_measured)  
}

load_marker_information <- function(samples,
                                    excel_file,
                                    sheet_name,
                                    recompute = FALSE){
  
  marker_info <- load_marker_excel(excel_file, sheet_name)
  
  markers_per_file <- get_markers_per_file(samples,
                                           result_rds = "RDS/markers_per_file.RDS",
                                           recompute)
  
  markers <- sort(unique(unlist(markers_per_file)))
  channels <- sapply(markers, function(m) return(NA))
  for(sample in names(markers_per_file)){
    channels[markers_per_file[[sample]]] <- names(markers_per_file[[sample]])
  }
  
  # Extra information about the markers
  marker_annotation <- data.frame(Marker = markers,
                                  Marker_short = gsub("^[^_]*_", "", markers),
                                  Channel = channels,
                                  CellType = marker_info[markers, "cell type"],
                                  MonocytePanel = marker_info[markers, "monocyte panel"],
                                  CyTOFPanel = marker_info[markers, "CyTOF panel"],
                                  MarkerType = marker_info[markers, "marker type"],
                                  stringsAsFactors = FALSE)
  rownames(marker_annotation) <- markers
  marker_annotation <- marker_annotation[order(marker_annotation$CyTOFPanel,
                                               marker_annotation$MarkerType,
                                               marker_annotation$CellType), ]
  
  # Which samples have which  markers
  marker_overview <- matrix(FALSE, 
                            nrow = nrow(samples),
                            ncol = nrow(marker_annotation),
                            dimnames = list(rownames(samples),
                                            rownames(marker_annotation)))
  for(sample in rownames(samples)) {
    marker_overview[sample, markers_per_file[[sample]]] <- TRUE
  }
  
  return(cbind(marker_annotation,
               t( marker_overview)[rownames(marker_annotation), ]))
}
