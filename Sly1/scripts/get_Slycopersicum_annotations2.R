library(biomartr)
library(tidyverse)
library("biomaRt")
library(org.At.tair.db)
library(GenomicFeatures)
library(tidyr)

get_Slycopersicum_annotations2 <- function(organism = "Solanum lycopersicum",
                                           attr = c("description",
                                                    "athaliana_eg_homolog_ensembl_gene",
                                                    "athaliana_eg_homolog_associated_gene_name",
                                                    "external_gene_name",
                                                    "athaliana_eg_homolog_perc_id")
) {
  all_slyc_genes <- GenomicFeatures::makeTxDbFromGFF("inputs/ITAG3.2_gene_models.gff", format = "gff3")
  all_slyc_genes <- genes(all_slyc_genes) %>% as.data.frame()
  
  all_slyc_genes_annotated <- biomartr::biomart(genes = all_slyc_genes$gene_id,
                                                mart = "plants_mart",                 
                                                dataset = "slycopersicum_eg_gene",           
                                                attributes = attr,
                                                filters = "ensembl_gene_id")
  
  all_slyc_genes_annotated$gene_id_versionless <- substr(all_slyc_genes_annotated$ensembl_gene_id, 1, 14)
  
  return(all_slyc_genes_annotated)
}
