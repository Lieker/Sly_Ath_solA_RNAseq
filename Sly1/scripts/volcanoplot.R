library(apeglm)
library(EnhancedVolcano)
source("Sly1/scripts/get_unfiltered_res_dds.R")
source("Sly1/scripts/get_DESeq_dds.R")
source("Sly1/scripts/coefficient_for_volcanoplot.R")

make_volcanoplot <- function(counts_csv_file = "Sly1/input/counts.csv",
                             xp_design_csv_file = "Sly1/input/xp_design.csv",
                             plantpart = "root",
                             ref_treatment = "no_solA",
                             treatment2 = "yes_solA",
                             log2FC_threshold = 0,
                             FCcutoff_volcano = 0,
                             padj_threshold = 0,
                             ttl = "",
                             method
) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       plantpart,
                       ref_treatment,
                       treatment2,
                       method)
  res <- get_unfiltered_res_dds(counts_csv_file,
                                xp_design_csv_file,
                                plantpart,
                                ref_treatment,
                                treatment2,
                                method)
  coef <- coefficient(counts_csv_file,
                      xp_design_csv_file,
                      plantpart,
                      ref_treatment,
                      treatment2,
                      method)
  
  shrunk <- lfcShrink(dds = dds,
                      res = res,
                      type = "apeglm",
                      coef = coef)
  v <- EnhancedVolcano(toptable = shrunk,
                       x = "log2FoldChange",
                       y = "padj",
                       pCutoff = padj_threshold,
                       FCcutoff = FCcutoff_volcano,
                       lab = rownames(shrunk),
                       title = ttl,
                       subtitle = "")
  return(v)
}

