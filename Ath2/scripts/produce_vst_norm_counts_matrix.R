suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("scripts/get_xp_design.R")

# this functions returns a matrix of scaled counts.
# samples in rows
# genes in columns
# row names = samples
produce_vst_norm_counts_matrix <- function(counts_csv_file = "inputs/counts.csv") {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>% 
    column_to_rownames("Geneid")
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  dds <- estimateSizeFactors(dds)
  dds_vst <- vst(dds)
  
  counts_vst <- assay(dds_vst, normalized = TRUE)
    
  counts_vst_t <-  t(counts_vst) # to have sample in rows
  counts_vst_t = as.data.frame(counts_vst_t)
  counts_vst_t <- counts_vst_t %>% rownames_to_column("sample")
  return(counts_vst_t)
}

