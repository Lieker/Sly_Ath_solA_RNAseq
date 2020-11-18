library(apeglm)
library(EnhancedVolcano)
source("scripts/DEG_function.R")

make_volcanoplot <- function(counts_csv_file = "raw_counts.csv",
                             xp_design_csv_file = "xp_design.csv",
                             timepoint = 2,
                             ref_treatment = "ethanol",
                             treatment2 = "millimolar_solanoeclepinA",
                             log2FC_threshold = 0,
                             padj_threshold = 0.05,
                             
) {
  list <- get_list_of_DEGs(counts_csv_file,
                           xp_design_csv_file,
                           timepoint,
                           ref_treatment,
                           treatment2,
                           log2FC_threshold,
                           padj_threshold)
  shrunk <- lfcShrink(dds = 
                        NOG NIET AF)
}


resultsNames(dds_2m)
res_2m_shr <- lfcShrink(dds = dds_2m,
                        res = res_2m,
                        type = "apeglm",
                        coef = "treatment_m_vs_e")

EnhancedVolcano(toptable = res_2m_shr, 
                x = "log2FoldChange",           
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                lab = rownames(res_2m_shr),
                title = 'DEGs t=2 1.5uM solA treatment',
                subtitle = 'LFC>1, padj <0.05')