suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))

# this functions returns a matrix of scaled counts.
# samples in rows
# genes in columns
# row names = samples
produce_scaled_counts_matrix <- function(count_csv_file = "raw_counts.csv", 
                                         xp_design_csv_file = "xp_design.csv") {
  counts <- read.csv(file = count_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>% 
    column_to_rownames("Geneid")
  
  xp_design <- read.csv(file = xp_design_csv_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  dds <- estimateSizeFactors(dds)
  dds_vst <- vst(dds)
  
  counts_vst <- assay(dds_vst, normalized = TRUE)
    
  counts_vst_t <-  t(counts_vst) # to have sample in rows
  counts_vst_t = as.data.frame(counts_vst_t) %>% 
    rownames_to_column("sample")
  
  return(counts_vst_t)
}