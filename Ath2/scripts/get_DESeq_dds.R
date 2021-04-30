library(DESeq2)
library(tidyr)
library(tibble)
source("Ath2/scripts/get_xp_design.R")

get_DESeq_dds <- function(counts_csv_file = "Ath2/input/counts.csv",
                          xp_design_csv_file = "Ath2/input/xp_design.csv",
                          trtm = c("a","b"),
                          ref_treatment = "a",
                          treatment2 = "b",
                          method = "treatment" #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                          ) {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
    column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file)
  
  if(method == "treatment"){
    xp_design <- xp_design[which(xp_design$treatment %in% trtm),] %>% 
      filter(treatment == ref_treatment | treatment == treatment2)
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  } else if (method == "solA"){
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ N + P + solA)
  } else if (method == "N:solA+solA"){
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ N + P + N:solA + solA)
  } else if (method == "N:solA"){
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ N + P + N:solA)
  } else if (method == "Ath12"){
    counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ experiment + treatment)
  } else{print("no valid method selected")}

  dds <- DESeq(dds)
  return(dds)
}
