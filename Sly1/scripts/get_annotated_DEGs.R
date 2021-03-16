library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(biomartr)


source("scripts/get_filtered_list_of_DEGs.R")
source("scripts/get_Slycopersicum_annotations.R")

get_annotated_DEGs <- function(counts_csv_file = "inputs/raw_counts.csv",
                               xp_design_csv_file = "inputs/xp_design.csv",
                               plantpart = "root",
                               ref_treatment = "no_solA",
                               treatment2 = "yes_solA",
                               log2FC_threshold = 0,
                               padj_threshold = 0.05,
                               name = "outputs/annotated_DEGslist.csv",
                               attr = c("description",
                                        "athaliana_eg_homolog_ensembl_gene",
                                        "athaliana_eg_homolog_associated_gene_name",
                                        "external_gene_name",
                                        "athaliana_eg_homolog_perc_id")
                               ) {
  res <- get_list_of_DEGs(counts_csv_file, 
                         xp_design_csv_file, 
                         plantpart, 
                         ref_treatment, 
                         treatment2, 
                         log2FC_threshold, 
                         padj_threshold) %>% rownames_to_column("gene")
   
  annotations <- get_Slycopersicum_annotations(attr = attr)
  res_annotated <- left_join(res, annotations, by = "gene")
  write_delim(x = res_annotated,
              file = name,
              delim = ";")
  
  return(res_annotated)
}
