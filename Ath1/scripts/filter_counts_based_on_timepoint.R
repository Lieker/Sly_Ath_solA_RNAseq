suppressPackageStartupMessages(library(tidyverse))
scaled_counts <- produce_scaled_counts_matrix()


filter_counts_based_on_timepoint <- function(counts = scaled_counts, 
                                        xp_design_csv_file = "Ath1/input/xp_design.csv", 
                                        tp = c(2,6,24)) {
  # read and filter xp design info for selected time
  xp_design = read.csv(file = xp_design_csv_file,
                       header = TRUE, 
                       check.names = FALSE, 
                       fileEncoding = "UTF-8-BOM") %>% dplyr::filter(time %in% tp)
  xp_design$time <- as.character(xp_design$time)
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  filtered_counts <- counts[counts$sample %in% xp_design$sample,] 
  
  row.names(filtered_counts) <- filtered_counts$sample
  filtered_counts$sample <- NULL
  
  return(filtered_counts)
}

