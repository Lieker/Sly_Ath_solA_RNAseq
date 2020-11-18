library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(org.At.tair.db)

source("scripts/DEG_function.R")
source("scripts/get_all_annotations.R")

get_annotated_DEGs <- function(counts_csv_file = "raw_counts.csv",
                               xp_design_csv_file = "xp_design.csv",
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
                               name = "annotated_DEGslist.csv") {
  
  res <- get_list_of_DEGs(counts_csv_file, 
                        xp_design_csv_file, 
                        timepoint, 
                        ref_treatment, 
                        treatment2, 
                        log2FC_threshold, 
                        padj_threshold) %>% rownames_to_column("genes")
  
  subset_annotated <- biomartr::biomart(genes = res$genes,
                                   mart       = "plants_mart",                 
                                   dataset    = "athaliana_eg_gene",           
                                   attributes = attr,        
                                   filters =  "ensembl_gene_id" )
  names(subset_annotated)[names(subset_annotated) == 'ensembl_gene_id'] <- 'genes'
  res_annotated <- left_join(res, subset_annotated, by = "genes")
  write_delim(x = res_annotated,
              file = name,
              delim = ";")
  
  return(res_annotated)
}

get_annotated_DEGs()
