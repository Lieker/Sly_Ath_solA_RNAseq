library(clusterProfiler)
library(org.At.tair.db)
source("Ath2/scripts/get_annotated_DEGs.R")
source("Ath2/scripts/get_Athaliana_annotations.R")

GO_annotation_list <- function(counts_csv_file = "Ath2/inputs/counts.csv",
                               xp_design_csv_file = "Ath2/inputs/xp_design.csv",
                               trtm = c("a","b"),
                               ref_treatment = "a",
                               treatment2 = "b",
                               method = "treatment", #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                               log2FC_threshold = 0,
                               padj_threshold = 0.05,
                               organism = "Arabidopsis thaliana",
                               attr = c("tair_symbol", 
                                        "uniprotswissprot",
                                        "entrezgene_id",
                                        "description",
                                        "external_gene_name",
                                        "external_gene_source"),
                               name = "Ath2/outputs/annotated_DEGslist.csv",
                               method2 = "DEG", #this parameter chooses if DEGs will be analysed or a comparison with LRT will be made
                               ont = "BP", #can be either BP, CC or MF
                               namego = "Ath2/outputs/GOannotationslist.csv",
                               plottype,
                               n = 30
                               ) {
  all_Ath_genes_annotated <- get_Athaliana_annotations(organism, attr, counts_csv_file)
  res_annotated <- get_annotated_DEGs(counts_csv_file, 
                                      xp_design_csv_file, 
                                      trtm, 
                                      ref_treatment, 
                                      treatment2,
                                      method,
                                      log2FC_threshold, 
                                      padj_threshold, 
                                      organism, 
                                      attr, 
                                      name,
                                      method2)
  go <- enrichGO(gene = res_annotated$entrezgene_id, 
                 universe = all_Ath_genes_annotated$entrezgene_id, 
                 OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                 keyType = "ENTREZID",
                 ont = ont,              # either "BP", "CC" or "MF",
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05,
                 readable = TRUE, 
                 pool = FALSE)
  write_delim(x = as.data.frame(go@result), 
              file = namego, 
              delim = ";")
  if(plottype == "bar"){
    d <- barplot(go, showCategory = n)
  } else if(plottype == "go"){
    d <- goplot(go, showCategory = n)
  }
  return(d)

}
