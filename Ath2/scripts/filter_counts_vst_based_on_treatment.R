suppressPackageStartupMessages(library(tidyverse))
source("scripts/get_xp_design.R")
source("scripts/produce_vst_norm_counts_matrix.R")

filter_counts_vst_based_on_treatment <- function(counts_csv_file = "inputs/counts.csv",
                                                 trtm = c("a","b","c","d","e","f","g")) {
  # filter xp_design info for selected time
  xp_design <- get_xp_design()
  xp_design <- xp_design[which(xp_design$treatment %in% trtm),]
  
  # filter counts based on xp_design info
  # place sample in row names for downstream compatibility with mypca() function
  counts_vst_t <- produce_vst_norm_counts_matrix(counts_csv_file)
  filtered_counts_vst_t <- counts_vst_t[counts_vst_t$sample %in% xp_design$sample,] 
  row.names(filtered_counts_vst_t) <- NULL
  filtered_counts_vst_t <- filtered_counts_vst_t %>% column_to_rownames("sample")
  
  
  return(filtered_counts_vst_t)
}
