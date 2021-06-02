library(DESeq2)
library(tidyr)
library(tibble)
source("Ath1/scripts/get_xp_design.R")

get_DESeq_dds <- function(counts_csv_file = "Ath1/input/counts.csv",
                          xp_design_csv_file = "Ath1/input/xp_design.csv",
                          trtm = c("ethanol","millimolar_solanoeclepinA"),
                          ref_treatment = "ethanol",
                          treatment2 = "millimolar_solanoeclepinA",
                          method = "time+treatment",
                          tp = c(2, 6, 24)
                          ) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
    column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file)
  
  if(method == "time+treatment"){
    xp_design <- xp_design[which(xp_design$treatment %in% trtm),] %>% 
      filter(treatment == ref_treatment | treatment == treatment2)
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ rna_isolation_batch + time + treatment)
  } else if(method == "treatment"){
    xp_design <- xp_design[which(xp_design$treatment %in% trtm),] %>%
      filter(treatment == ref_treatment | treatment == treatment2)
    xp_design <- xp_design[which(xp_design$time %in% tp),]
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ rna_isolation_batch + treatment)
  }
  
  dds <- DESeq(dds)
  return(dds)
}
