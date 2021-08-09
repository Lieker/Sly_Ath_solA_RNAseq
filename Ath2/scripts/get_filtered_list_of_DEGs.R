library(DESeq2)
library(dplyr)
library(tidyr)
source("Ath2/scripts/get_unfiltered_res_dds.R")

get_list_of_DEGs <- function(d= dds,
                             r = res,
                             method = "treatment", #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                             log2FC_threshold = 0,
                             padj_threshold = 0.05) {
  
  if(method == "treatment"){
    res_shr <- lfcShrink(dds = d,
                       res = r,
                       type = "apeglm",
                       coef = resultsNames(d)[2])
    
  } else if(method == "solA") {
    res_shr <- lfcShrink(dds = d,
                         res = r,
                         type = "apeglm",
                         coef = resultsNames(d)[4])
  }
  
  res_shr <- res_shr %>% as.data.frame() %>% filter(padj < padj_threshold) %>% filter(log2FoldChange > log2FC_threshold | -log2FoldChange > log2FC_threshold)
  
  return(res_shr)
}
