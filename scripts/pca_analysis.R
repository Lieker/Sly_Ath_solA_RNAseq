suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
source("scripts/pca_function.R")


#####################
# Import data
#####################
counts <- read.csv("raw_counts.csv", header = TRUE, stringsAsFactors = FALSE) %>% 
  column_to_rownames("Geneid") %>% 
  t() # so that samples are in rows
  
xp_design <- read.csv("xp_design.csv", header = TRUE, stringsAsFactors = FALSE)



#####
# PCA 
#####
pca_res <- mypca(counts, center = TRUE, scale = TRUE)

scree_plot <- ggplot(pca_res$explained_var, aes(x = PC, y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('explained variance per component') + 
  geom_bar(stat = "identity")
scree_plot


