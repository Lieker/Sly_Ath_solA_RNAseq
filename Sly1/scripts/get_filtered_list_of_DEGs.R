library(DESeq2)
library(dplyr)
library(tidyr)
source("scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "inputs/raw_counts.csv",
                             xp_design_csv_file = "inputs/xp_design.csv",
                             plantpart = "root",
                             ref_treatment = "no_solA",
                             treatment2 = "yes_solA",
                             log2FC_threshold = 1,
                             padj_threshold = 0.05) {
  res_filtered <- get_unfiltered_res_dds(counts_csv_file, 
                                         xp_design_csv_file,
                                         plantpart,
                                         ref_treatment,
                                         treatment2) 
  res_filtered <- as.data.frame(res_filtered) %>%
    filter(padj < padj_threshold) %>% 
    filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_filtered)
}
r <- get_list_of_DEGs()
