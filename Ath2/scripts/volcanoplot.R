library(apeglm)
library(EnhancedVolcano)
source("scripts/get_unfiltered_res_dds.R")
source("scripts/get_DESeq_dds.R")
source("scripts/coefficient_for_volcanoplot.R")

make_volcanoplot <- function(counts_csv_file = "inputs/counts.csv",
                             xp_design_csv_file = "inputs/xp_design.csv",
                             trtm = c("a","b"),
                             ref_treatment = "a",
                             treatment2 = "b",
                             log2FC_threshold = 0,
                             FCcutoff_volcano = 0,
                             padj_threshold = 0
) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2)
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                trtm,
                                ref_treatment,
                                treatment2)
  coef <- coefficient(counts_csv_file,
                      xp_design_csv_file,
                      trtm,
                      ref_treatment,
                      treatment2)
  
  shrunk <- lfcShrink(dds = dds,
                      res = res,
                      type = "apeglm",
                      coef = coef)
  v <- EnhancedVolcano(toptable = shrunk,
                       x = "log2FoldChange",
                       y = "padj",
                       pCutoff = padj_threshold,
                       FCcutoff = FCcutoff_volcano,
                       lab = rownames(shrunk) )
  return(v)
}
