library(DESeq2)
library(tidyr)
library(tibble)
source("Sly1/scripts/get_xp_design.R")

get_DESeq_dds <- function(counts_csv_file = "Sly1/input/counts.csv",
                          xp_design_csv_file = "Sly1/input/xp_design.csv",
                          plantpart = "root",
                          ref_treatment = "no_solA",
                          treatment2 = "yes_solA",
                          method) 
  {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE,
                     fileEncoding = "UTF-8-BOM") %>% column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file) %>%
    filter(compartment == plantpart) %>% 
    filter(treatment == ref_treatment | treatment == treatment2)
  
  counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
  
  if(method == "treatment") {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = xp_design,
                                design = ~ treatment)
  } else if(method == "compartment+treatment") {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = xp_design,
                                  design = ~ compartment + treatment)
  }
  
  dds <- DESeq(dds)
  
  return(dds)
}
