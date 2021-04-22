library(apeglm)
library(EnhancedVolcano)
source("Ath2/scripts/get_unfiltered_res_dds.R")
source("Ath2/scripts/get_DESeq_dds.R")
source("Ath2/scripts/coefficient_for_volcanoplot.R")

make_volcanoplot <- function(counts_csv_file = "Ath2/inputs/counts.csv",
                             xp_design_csv_file = "Ath2/inputs/xp_design.csv",
                             trtm = c("a","b"),
                             ref_treatment = "a",
                             treatment2 = "b",
                             method = "treatment", #this parameter chooses which formula design will be chosen: ~treatment or ~N + P + solA
                             log2FC_threshold = 0,
                             FCcutoff_volcano = 0,
                             padj_threshold = 0,
                             xsize = 2,
                             ttl = ""
) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method)
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                trtm,
                                ref_treatment,
                                treatment2,
                                method)
  coeff <- coefficient(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method)
  
  shrunk <- lfcShrink(dds = dds,
                      res = res,
                      type = "apeglm",
                      coef = coeff)
  v <- EnhancedVolcano(toptable = shrunk,
                       x = "log2FoldChange",
                       y = "padj",
                       title = ttl,
                       subtitle = "",
                       pCutoff = padj_threshold,
                       FCcutoff = FCcutoff_volcano,
                       lab = rownames(shrunk),
                       xlim = c(-(xsize), xsize))
  return(v)
}
