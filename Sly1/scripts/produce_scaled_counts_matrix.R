suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))

# this functions returns a matrix of scaled counts.
# samples in rows
# genes in columns
# row names = samples
produce_scaled_counts_matrix <- function(count_csv_file = "inputs/raw_counts.csv", 
                                         xp_design_csv_file = "inputs/xp_design.csv") {
  counts <- read.csv(file = count_csv_file, 
                     header = TRUE, 
                     stringsAsFactors = FALSE) %>% 
    column_to_rownames("Geneid")
  
  xp_design <- read.csv(file = xp_design_csv_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ 1)
  dds <- estimateSizeFactors(dds)
  dds = estimateDispersions(object = dds, 
                            fitType = "parametric", 
                            quiet = TRUE)
  vsd = varianceStabilizingTransformation(object = dds, 
                                          blind = TRUE,           # do not take the design formula into account. 
                                          fitType = "parametric") # best practice for sample-level QC
                                          
  
  variance_stabilised_counts <- assay(vsd)
  
  variance_stabilised_counts = t(variance_stabilised_counts) # to have sample in rows
  scaled_counts = as.data.frame(variance_stabilised_counts) %>% 
    rownames_to_column("sample")
  
  return(scaled_counts)
}
