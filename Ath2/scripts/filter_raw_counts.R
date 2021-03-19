library(tidyr)
library(dplyr)
library(tidyverse)

filter_raw_counts <- function(counts_csv_file = 'inputs/counts.csv', 
                              trtm = c("a","b","c","d","e","f","g")) {
  # read and filter xp design info for selected time
  xp_design <- get_xp_design()
  xp_design_filtered <- xp_design[which(xp_design$treatment %in% trtm),]
  
  counts <- counts <- read.csv(file = counts_csv_file, 
                               header = TRUE, 
                               stringsAsFactors = FALSE) %>% column_to_rownames("Geneid") %>%  t() %>% as.data.frame() %>% rownames_to_column("sample")
  
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  
  filtered_counts <- counts[counts$sample %in% xp_design_filtered$sample,] %>% column_to_rownames("sample") %>% t()
  
  
  return(filtered_counts)
}
