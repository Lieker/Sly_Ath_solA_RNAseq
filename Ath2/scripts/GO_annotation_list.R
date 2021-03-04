library(clusterProfiler)
library(org.At.tair.db)
source("get_annotated_DEGs.R")
source("scripts/get_Athaliana_annotations")

GO_annotation_list <- function(counts_csv_file = "inputs/raw_counts.csv",
                               xp_design_csv_file = "inputs/xp_design.csv",
                               timepoint = 2,
                               ref_treatment = "ethanol",
                               treatment2 = "millimolar_solanoeclepinA",
                               log2FC_threshold = 0,
                               padj_threshold = 0.05,
                               organism = "Arabidopsis thaliana",
                               attr = c("tair_symbol", 
                                        "uniprotswissprot",
                                        "entrezgene_id",
                                        "description",
                                        "external_gene_name",
                                        "external_gene_source"),
                               name = "annotated_DEGslist.csv",
                               ont = "BP", #can be either BP, CC or MF
                               namego = "outputs/GOannotationslist.csv"
                               ) {
  all_Ath_genes_annotated <- get_Athaliana_annotations(organism, attr, counts_csv_file)
  res_annotated <- get_annotated_DEGs(counts_csv_file, 
                                      xp_design_csv_file, 
                                      timepoint, 
                                      ref_treatment, 
                                      treatment2, 
                                      log2FC_threshold, 
                                      padj_threshold, 
                                      organism, 
                                      attr, 
                                      name)
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
  d <- dotplot(go)
  return(d)
}
