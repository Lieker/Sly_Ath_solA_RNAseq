library(DESeq2)
library(tidyr)
library(tibble)
source("scripts/get_xp_design.R")

get_DESeq_dds <- function(counts_csv_file = "inputs/counts.csv",
                          xp_design_csv_file = "inputs/xp_design.csv",
                          trtm = c("a","b"),
                          ref_treatment = "a",
                          treatment2 = "b") {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
    column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file)
  xp_design <- xp_design[which(xp_design$treatment %in% trtm),] %>% 
    filter(treatment == ref_treatment | treatment == treatment2)
  
  counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ N + treatment)
  dds <- DESeq(dds)
  return(dds)
}
