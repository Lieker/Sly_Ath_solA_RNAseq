compare_wald_vs_LRT <- function(counts_csv_file = "Sly1/input/counts.csv",
                                xp_design_csv_file = "Sly1/input/xp_design.csv",
                                compartment = "root",
                                ref_treatment = "no_solA",
                                treatment2 = "yes_solA",
                                method, #this parameter chooses which formula design will be used
                                padj_cutoff = 0.05
)
{
  source("Ath2/scripts/get_DESeq_dds.R")
  dds <- get_DESeq_dds(counts_csv_file = counts_csv_file,
                       xp_design_csv_file = xp_design_csv_file,
                       trtm = trtm,
                       ref_treatment = ref_treatment,
                       treatment2 = treatment2,
                       method = method)
  if(method == "treatment"){
    dds_lrt <- DESeq(dds, test="LRT", reduced = ~1)
  } else if(method == "solA" | method == "N:solA"){
    dds_lrt <- DESeq(dds, test="LRT", reduced = ~ N + P)
  } else if(method == "N:solA+solA") {
    dds_lrt <- DESeq(dds, test="LRT", reduced = ~ N + P + N:solA)
  }
  
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

  overlap <- sig_res_wald[sig_res_wald$gene %in% sig_res_LRT$gene,]
  
  print(paste0("LRT gives ", nrow(sig_res_LRT), " DEGs, while Wald test gives ", nrow(sig_res_wald), " DEGs."))
  
  return(overlap)
}
