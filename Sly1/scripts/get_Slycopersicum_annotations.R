library(biomartr)
library(tidyverse)
library("biomaRt")
library(org.At.tair.db)
library(stringr)

get_Slycopersicum_annotations <- function(attr = c("description",
                                                   "athaliana_eg_homolog_ensembl_gene",
                                                   "athaliana_eg_homolog_associated_gene_name",
                                                   "external_gene_name",
                                                   "athaliana_eg_homolog_perc_id") 
) {
  descriptions <- read.csv("inputs/ITAG4.1_descriptions.csv", stringsAsFactors = FALSE)
  descriptions$gene <- substr(descriptions$gene, 1, 16)
  
  go_terms <- read.csv("inputs/ITAG4.1_goterms.csv", stringsAsFactors = FALSE)
  go_terms$gene <- substr(go_terms$gene, 1, 16)
  
  annotations <- left_join(descriptions, go_terms, by = "gene")
  annotations$gene_id_versionless <- substr(annotations$gene, 1, 14)
  
  source("scripts/get_Slycopersicum_annotations2.R")
  
  ensembl_annotations <- get_Slycopersicum_annotations2(attr = attr)
  
  annotations <- left_join(annotations, ensembl_annotations, by = "gene_id_versionless", suffix = c("_ITAG4.0", "_ITAG3.0"))
  
  return(annotations)
}
  



