suppressPackageStartupMessages(library(tidyverse))


filter_counts_based_on_treatment <- function(counts = scaled_counts, 
                                        xp_design_csv_file = "inputs/xp_design.csv", 
                                        tr = c("a","b","c","d","e","f","g")) {
  # read and filter xp design info for selected time
  xp_design = read.csv(file = xp_design_csv_file,
                       header = TRUE, 
                       check.names = FALSE, 
                       fileEncoding = "UTF-8-BOM") %>% dplyr::filter(treatment %in% tr)
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  filtered_counts <- counts[counts$sample %in% xp_design$sample,] 
  
  row.names(filtered_counts) <- filtered_counts$sample
  filtered_counts$sample <- NULL
  
  return(filtered_counts)
}

