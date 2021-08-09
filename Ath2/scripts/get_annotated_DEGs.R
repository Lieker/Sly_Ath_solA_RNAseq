library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(org.At.tair.db)

source("Ath2/scripts/compare_Wald_vs_LRT.R")
source("Ath2/scripts/get_filtered_list_of_DEGs.R")
source("Ath2/scripts/get_Athaliana_annotations.R")

get_annotated_DEGs <- function(d=dds,
                               r=res,
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
                               name = "Ath2/output/annotated_DEGslist.csv",
                               method2 = "DEG" #this parameter chooses if DEGs acc. to Wald will be analysed or a comparison with LRT will be made
                               ) {
  if (method2 == "DEG") {
    res <- get_list_of_DEGs(d,r,
                            method,
                            log2FC_threshold, 
                            padj_threshold) %>% rownames_to_column("gene")
  } else if (method2 == "LRTcompare") {
    res <- compare_wald_vs_LRT(trtm = trtm, ref_treatment = ref_treatment, treatment2 = treatment2)
  } else { 
    print("no valid method2; choose DEG or LRTcompare") }
  
  genes <- res$gene
  names(res)[names(res) == "gene"] <- "transcript"
  genes <- substr(genes, start = 1, stop = 9)
  res$gene <- genes

  subset_annotated <- biomartr::biomart(genes = res$gene,
                                        mart = "plants_mart",                 
                                        dataset = "athaliana_eg_gene",           
                                        attributes = attr,        
                                        filters =  "ensembl_gene_id" )
  names(subset_annotated)[names(subset_annotated) == 'ensembl_gene_id'] <- 'gene'
  
  res_annotated <- left_join(res, subset_annotated, by = "gene")
  uniquepaste <- function(x) { paste(unique(x), sep = ',', collapse = ",") }
  res_annotated <- aggregate(res_annotated, by = list(res_annotated$gene, res_annotated$baseMean, res_annotated$log2FoldChange), FUN = uniquepaste)
  write_delim(res_annotated, path = name, delim = ",", col_names = TRUE)
  
  return(res_annotated)
}

