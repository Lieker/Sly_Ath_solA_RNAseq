library(biomartr)
library(tidyverse)
library("biomaRt")
library(org.At.tair.db)

get_Athaliana_annotations <- function(organism = "Arabidopsis thaliana",
                                      attr = c("tair_symbol", 
                                      "uniprotswissprot",
                                      "entrezgene_id",
                                      "description",
                                      "external_gene_name",
                                      "external_gene_source"),
                                      counts_csv_file = "inputs/counts.csv"
) {biomartr::organismBM(organism = organism)
  
  arabido_attributes = biomartr::organismAttributes(organism) %>% 
  filter(dataset == "athaliana_eg_gene")
  
  all_arabidopsis_genes <- read.csv(counts_csv_file, header = TRUE, stringsAsFactors = FALSE)
  all_arabidopsis_genes <- all_arabidopsis_genes$Geneid
  all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                       mart       = "plants_mart",                 
                                                       dataset    = "athaliana_eg_gene",           
                                                       attributes = attr,        
                                                       filters =  "ensembl_gene_id" )  
  all_arabidopsis_genes_annotated$entrezgene_id = as.character(all_arabidopsis_genes_annotated$entrezgene_id) 
  return(all_arabidopsis_genes_annotated)
}
  



