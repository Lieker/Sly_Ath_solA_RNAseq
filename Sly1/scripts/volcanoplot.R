library(apeglm)
library(EnhancedVolcano)
source("Sly1/scripts/get_unfiltered_res_dds.R")
source("Sly1/scripts/get_DESeq_dds.R")

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
  
  
  shrunk <- lfcShrink(dds = dds,
                      res = res,
                      type = "apeglm",
                      coef = resultsNames(dds)[2])
  v <- EnhancedVolcano(toptable = shrunk,
                       x = "log2FoldChange",
                       y = "padj",
                       pCutoff = padj_threshold,
                       FCcutoff = FCcutoff_volcano,
                       lab = rownames(shrunk),
                       selectLab = "",
                       title = ttl,
                       subtitle = "",
                       xlim = c(-7, 7),
                       ylim = c(0, -log10(10e-4)),
                       colGradient = c("red3", "royalblue"))
  return(v)
}

