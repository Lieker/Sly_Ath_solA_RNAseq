suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("scripts/pca_function.R")


#############################
# Scale counts the DESeq2 way
#############################
counts <- read.csv("raw_counts.csv", header = TRUE, stringsAsFactors = FALSE) %>% 
  column_to_rownames("Geneid")

xp_design <- read.csv("xp_design.csv", header = TRUE, stringsAsFactors = FALSE)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = xp_design, design = ~ 1)
dds <- estimateSizeFactors(dds)
scaled_counts = counts(dds, normalized = TRUE)

scaled_counts = t(scaled_counts) # to have sample in rows

#####
# PCA 
#####
pca_res <- mypca(counts, center = TRUE, scale = TRUE)

scree_plot <- ggplot(pca_res$explained_var, aes(x = PC, y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component') + 
  geom_bar(stat = "identity")
scree_plot


