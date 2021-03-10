library(biomartr)
library(tidyverse)
library("biomaRt")
library(org.At.tair.db)

get_Slycopersicum_annotations <- function(organism = "Solanum lycopersicum",
                                      attr = c("description",
                                               "athaliana_eg_homolog_ensembl_gene",
                                               "athaliana_eg_homolog_associated_gene_name",
                                               "external_gene_name"),
                                      counts_csv_file = "raw_counts.csv"
) {biomartr::organismBM(organism = organism)
  
  attributes = biomartr::organismAttributes(organism) %>% 
  filter(dataset == "slycopersicum_eg_gene")
  
  all_slyc_genes <- read.csv(counts_csv_file, header = TRUE, stringsAsFactors = FALSE)
  all_slyc_genes <- all_slyc_genes$Geneid
  all_slyc_genes_annotated <- biomartr::biomart(genes = all_slyc_genes,
                                                       mart       = "plants_mart",                 
                                                       dataset    = "slycopersicum_eg_gene",           
                                                       attributes = attr,        
                                                       filters =  "ensembl_gene_id" )  
  all_slyc_genes_annotated$entrezgene_id = as.character(all_slyc_genes_annotated$entrezgene_id) 
  return(all_slyc_genes_annotated)
}
  



