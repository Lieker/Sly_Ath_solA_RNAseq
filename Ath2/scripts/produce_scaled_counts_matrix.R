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
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  dds = estimateDispersions(object = dds, 
                            fitType = "parametric", 
                            quiet = TRUE)
  vsd = varianceStabilizingTransformation(object = dds, 
                                          blind = TRUE,           # do not take the design formula into account. 
                                                                  # best practice for sample-level QC
                                          fitType = "parametric") # extract the matrix of variance stabilised counts
  variance_stabilised_counts <- assay(vsd)
  t_variance_stabilised_counts <- t(variance_stabilised_counts)
  t_variance_stabilised_counts <- as.data.frame(t_variance_stabilised_counts)
  
  return(t_variance_stabilised_counts)
}
