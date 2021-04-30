library(DESeq2)
library(dplyr)
library(tidyr)
source("Sly1/scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "Sly1/input/counts.csv",
                             xp_design_csv_file = "Sly1/input/xp_design.csv",
                             plantpart = "root",
                             ref_treatment = "no_solA",
                             treatment2 = "yes_solA",
                             log2FC_threshold = 1,
                             padj_threshold = 0.05,
                             method) {
  res_filtered <- get_unfiltered_res_dds(counts_csv_file, 
                                         xp_design_csv_file,
                                         plantpart,
                                         ref_treatment,
                                         treatment2,
                                         method) 
  res_filtered <- as.data.frame(res_filtered) %>%
    filter(padj < padj_threshold) %>% 
    filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_filtered)
}
