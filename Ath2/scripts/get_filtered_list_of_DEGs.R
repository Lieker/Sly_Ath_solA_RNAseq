library(DESeq2)
library(dplyr)
library(tidyr)
source("scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "inputs/counts.csv",
                             xp_design_csv_file = "inputs/xp_design.csv",
                             trtm = c("a","b"),
                             ref_treatment = "a",
                             treatment2 = "b",
                             log2FC_threshold = 1,
                             padj_threshold = 0.05) {
  res_filtered <- get_unfiltered_res_dds(counts_csv_file, 
                                         xp_design_csv_file,
                                         trtm,
                                         ref_treatment,
                                         treatment2) %>% as.data.frame() %>%
    filter(padj < padj_threshold) %>% 
    filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_filtered)
}