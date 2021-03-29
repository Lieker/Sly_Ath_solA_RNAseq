compare_wald_vs_LRT <- function(counts_csv_file = "inputs/counts.csv",
                                xp_design_csv_file = "inputs/xp_design.csv",
                                trtm = c("a","b"),
                                ref_treatment = "a",
                                treatment2 = "b",
                                padj_cutoff = 0.05
)
{
  source("scripts/get_DESeq_dds.R")
  dds <- get_DESeq_dds(counts_csv_file = counts_csv_file,
                       xp_design_csv_file = xp_design_csv_file,
                       trtm = trtm,
                       ref_treatment = ref_treatment,
                       treatment2 = treatment2)
  
  dds_lrt <- DESeq(dds, test="LRT", reduced = ~1)
  res_LRT <- results(dds_lrt)
  sig_res_LRT <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < padj_cutoff)
  sigLRT_genes <- sig_res_LRT %>% 
    pull(gene)
  
  res_wald <- results(dds)
  sig_res_wald <- res_wald %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>% 
    filter(padj < padj_cutoff)
  sigwald_genes <- sig_res_wald %>% 
    pull(gene)
  
  diff <- sig_res_wald[!(sig_res_wald$gene %in% sig_res_LRT$gene),]
<<<<<<< HEAD
  overlap <- sig_res_wald[sig_res_wald$gene %in% sig_res_LRT$gene,]
  
  print(paste0("LRT gives ", nrow(sig_res_LRT), " DEGs, while Wald test gives ", nrow(sig_res_wald), " DEGs."))
  
  return(overlap)
=======
  
  print(paste0("LRT gives ", nrow(sig_res_LRT), " DEGs, while Wald test gives ", nrow(sig_res_wald), " DEGs."))
  
  return(diff)
>>>>>>> 607b218894ec322033499f9291fc5c9e8c900ee3
  
}
