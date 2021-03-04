library(DESeq2)
library(dplyr)
library(tidyr)
source("scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "raw_counts.csv",
                             xp_design_csv_file = "xp_design.csv",
                             timepoint = 2,
                             ref_treatment = "ethanol",
                             treatment2 = "millimolar_solanoeclepinA",
                             log2FC_threshold = 1,
                             padj_threshold = 0.05) {
  res_filtered <- get_unfiltered_res_dds(counts_csv_file, 
                                         xp_design_csv_file,
                                         timepoint,
                                         ref_treatment,
                                         treatment2) %>% as.data.frame() %>%
    filter(padj < padj_threshold) %>% 
    filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_filtered)
}
r <- get_list_of_DEGs()
