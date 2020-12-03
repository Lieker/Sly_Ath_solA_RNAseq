library(DESeq2)
library(tidyr)
library(tibble)
source("scripts/get_xp_design.R")

get_DESeq_dds <- function(counts_csv_file = "raw_counts.csv",
                          xp_design_csv_file = "xp_design.csv",
                          timepoint = 2,
                          ref_treatment = "ethanol",
                          treatment2 = "millimolar_solanoeclepinA") {
  counts <- read.csv(file = counts_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>%
    column_to_rownames("Geneid")
  counts <- as.data.frame(t(counts)) %>% rownames_to_column("sample")
  
  xp_design <- get_xp_design(xp_design_csv_file) %>%
    filter(time == timepoint) %>% 
    filter(treatment == ref_treatment | treatment == treatment2)
  
  counts <- counts %>% filter(sample %in% xp_design$sample) %>% column_to_rownames("sample") %>% t()
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  dds <- DESeq(dds)
  return(dds)
}
