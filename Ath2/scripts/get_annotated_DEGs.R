library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(org.At.tair.db)

source("scripts/compare_Wald_vs_LRT.R")
source("scripts/get_filtered_list_of_DEGs.R")
source("scripts/get_Athaliana_annotations.R")

get_annotated_DEGs <- function(counts_csv_file = "inputs/counts.csv",
                               xp_design_csv_file = "inputs/xp_design.csv",
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
                               name = "outputs/annotated_DEGslist.csv",
                               method2 = "DEG" #this parameter chooses if DEGs will be analysed or a comparison with LRT will be made
                               ) {
  if (method2 == "DEG") {
    res <- get_list_of_DEGs(counts_csv_file, 
                            xp_design_csv_file, 
                            trtm, 
                            ref_treatment, 
                            treatment2,
                            method,
                            log2FC_threshold, 
                            padj_threshold) %>% rownames_to_column("gene")
  } else if (method2 == "LRTcompare") {
    res <- compare_wald_vs_LRT(trtm = trtm, ref_treatment = ref_treatment, treatment2 = treatment2)
  } else { 
    print("no valid method") }
  
  
  subset_annotated <- biomartr::biomart(genes = res$gene,
                                        mart = "plants_mart",                 
                                        dataset = "athaliana_eg_gene",           
                                        attributes = attr,        
                                        filters =  "ensembl_gene_id" )
  names(subset_annotated)[names(subset_annotated) == 'ensembl_gene_id'] <- 'gene'
  res_annotated <- left_join(res, subset_annotated, by = "gene")
  write_delim(x = res_annotated,
              file = name,
              delim = ";")
  
  return(res_annotated)
}
