library(DESeq2)
library(dplyr)
library(tidyr)
library(tidyverse)
source("scripts/get_xp_design.R")

get_list_of_DEGs <- function(counts_csv_file = "inputs/raw_counts.csv",
                             xp_design_csv_file = "inputs/xp_design.csv",
                             timepoint = 2,
                             ref_treatment = "ethanol",
                             treatment2 = "millimolar_solanoeclepinA",
                             log2FC_threshold = 0,
                             padj_threshold = 0.05) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
    column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file) %>% 
    dplyr::filter(time == timepoint) %>% 
    dplyr::filter(treatment == ref_treatment | treatment == treatment2)
  
  counts <- counts %>% dplyr::filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = xp_design, 
                                design = ~ rna_isolation_batch + treatment)
  dds$treatment <- relevel(dds$treatment, ref = ref_treatment)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds)) %>% 
    dplyr::filter(padj < padj_threshold) %>% 
    dplyr::filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res)
}


