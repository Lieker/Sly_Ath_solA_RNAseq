suppressPackageStartupMessages(library(tidyverse))


filter_counts_based_on_compartment <- function(counts = scaled_counts, 
                                               xp_design_csv_file = "inputs/xp_design.csv", 
                                               comprtmt = "root") {
  # read and filter xp design info for selected compartment
  xp_design = read.csv(file = xp_design_csv_file,
                       header = TRUE,
                       check.names = FALSE,
                       fileEncoding = "UTF-8-BOM") %>% 
    dplyr::filter(organ == comprtmt) %>% as.data.frame()
  
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  counts <- as.data.frame(counts)
  filtered_counts <- counts[counts$sample == xp_design$sample,] %>% as.data.frame
  
  
  row.names(filtered_counts) <- filtered_counts$sample
  filtered_counts$sample <- NULL
  
  return(filtered_counts)
}

