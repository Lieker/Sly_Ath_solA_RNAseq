suppressPackageStartupMessages(library(tidyverse))

source("scripts/produce_scaled_counts_matrix.R")

filter_counts_based_on_treatment <- function(counts = scaled_counts, 
                                             xp_design_csv_file = "inputs/xp_design.csv", 
                                             trtm = c("a","b","c","d","e","f","g")) {
  # read and filter xp design info for selected time
  xp_design <- get_xp_design()
  xp_design_filtered <- xp_design[which(xp_design$treatment %in% trtm),]
  
  counts <- produce_scaled_counts_matrix(count_csv_file = "inputs/counts.csv", 
                                         xp_design_csv_file = "inputs/xp_design.csv")
  
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  
  filtered_counts <- counts[counts$sample %in% xp_design_filtered$sample,] 
  
  row.names(filtered_counts) <- filtered_counts$sample
  filtered_counts$sample <- NULL
  
  return(filtered_counts)
}

