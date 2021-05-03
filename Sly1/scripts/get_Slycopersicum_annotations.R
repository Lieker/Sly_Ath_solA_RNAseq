get_Slycopersicum_annotations <- function(attr = c("description",
                                                   "athaliana_eg_homolog_ensembl_gene",
                                                   "athaliana_eg_homolog_associated_gene_name",
                                                   "external_gene_name",
                                                   "athaliana_eg_homolog_perc_id") 
) {
  get_Slycopersicum_annotations2 <- function(organism = "Solanum lycopersicum",
                                             attr = c("description",
                                                      "athaliana_eg_homolog_ensembl_gene",
                                                      "athaliana_eg_homolog_associated_gene_name",
                                                      "external_gene_name",
                                                      "athaliana_eg_homolog_perc_id")
  ) {
    all_slyc_genes <- GenomicFeatures::makeTxDbFromGFF("Sly1/input/ITAG3.2_gene_models.gff", format = "gff3")
    all_slyc_genes <- genes(all_slyc_genes) %>% as.data.frame()
    
    all_slyc_genes_annotated <- biomartr::biomart(genes = all_slyc_genes$gene_id,
                                                  mart = "plants_mart",                 
                                                  dataset = "slycopersicum_eg_gene",           
                                                  attributes = attr,
                                                  filters = "ensembl_gene_id")
    
    all_slyc_genes_annotated$gene_id_versionless <- substr(all_slyc_genes_annotated$ensembl_gene_id, 1, 14)
    
    return(all_slyc_genes_annotated)
  }
  
  descriptions <- read.csv("Sly1/input/ITAG4.1_descriptions.csv", stringsAsFactors = FALSE)
  descriptions$gene <- substr(descriptions$gene, 1, 16)
  
  go_terms <- read.csv("Sly1/input/ITAG4.1_goterms.csv", stringsAsFactors = FALSE)
  go_terms$gene <- substr(go_terms$gene, 1, 16)
  
  annotations <- left_join(descriptions, go_terms, by = "gene")
  annotations$gene_id_versionless <- substr(annotations$gene, 1, 14)
  
  ensembl_annotations <- get_Slycopersicum_annotations2(attr = attr)
  
  annotations <- left_join(annotations, ensembl_annotations, by = "gene_id_versionless", suffix = c("_ITAG4.0", "_ITAG3.0"))
  
  return(annotations)
}
  



