suppressPackageStartupMessages(library(tidyverse))

filter_counts_based_on_treatment <- function(counts = t_variance_stabilised_counts,
                                             xp_design_csv_file = "inputs/xp_design.csv", 
                                             trtm = c("a","b","c","d","e","f","g")) {
  # read and filter xp design info for selected time
  source("scripts/get_xp_design.R")
  xp_design <- get_xp_design()
  xp_design_filtered <- xp_design[which(xp_design$treatment %in% trtm),]
  source("scripts/produce_scaled_counts_matrix.R")
  counts <- produce_scaled_counts_matrix()
  counts <- as.data.frame(counts) %>% rownames_to_column("sample")
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  
  filtered_counts <- counts[counts$sample %in% xp_design_filtered$sample,] %>% column_to_rownames("sample")
  
  return(filtered_counts)
}

