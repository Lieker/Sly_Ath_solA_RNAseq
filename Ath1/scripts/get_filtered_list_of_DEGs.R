library(DESeq2)
library(dplyr)
library(tidyr)
source("Ath1/scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(counts_csv_file = "Ath1/input/counts.csv",
                             xp_design_csv_file = "Ath1/input/xp_design.csv",
                             trtm = c("ethanol","millimolar_solanoeclepinA"),
                             ref_treatment = "ethanol",
                             treatment2 = "millimolar_solanoeclepinA",
                             method = "time+treatment",
                             tp = c(2, 6, 24), 
                             log2FC_threshold = 0,
                             padj_threshold = 0.05) {
  res_filtered <- get_unfiltered_res_dds(counts_csv_file,
                                         xp_design_csv_file,
                                         trtm,
                                         ref_treatment,
                                         treatment2,
                                         method,
                                         tp) %>% as.data.frame() %>% dplyr::filter(padj < padj_threshold) %>% dplyr::filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_filtered)
}
