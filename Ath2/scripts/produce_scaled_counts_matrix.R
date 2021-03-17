suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))

# this functions returns a matrix of scaled counts.
# samples in rows
# genes in columns
# row names = samples
produce_scaled_counts_matrix <- function(count_csv_file = "inputs/counts.csv", 
                                         xp_design_csv_file = "inputs/xp_design.csv") {
  
  counts <- read.csv(file = count_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>% column_to_rownames("Geneid")
  
  xp_design <- read.csv(file = xp_design_csv_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  xp_design$treatment <- as.factor(xp_design$treatment)
  

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ treatment)
  dds <- estimateSizeFactors(dds)
  scaled_counts = counts(dds, normalized = TRUE)
  
  scaled_counts = t(scaled_counts) # to have sample in rows
  scaled_counts = as.data.frame(scaled_counts) %>% 
    rownames_to_column("sample")
  
  return(scaled_counts)
}
