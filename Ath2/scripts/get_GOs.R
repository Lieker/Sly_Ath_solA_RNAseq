library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(org.At.tair.db)

source("Ath2/scripts/compare_Wald_vs_LRT.R")
source("Ath2/scripts/get_filtered_list_of_DEGs.R")
source("Ath2/scripts/get_Athaliana_annotations.R")

get_GOs <- function(d=dds_all,
                    rlist,
                    log2FC_threshold = 0,
                    padj_threshold = 0.05,
                    organism = "Arabidopsis thaliana",
                    attr = c("tair_symbol",
                             "entrezgene_id",
                             "go_id",
                             "name_1006",
                             "definition_1006",
                             "go_linkage_type",
                             "namespace_1003")
                    ) {
  r1 <- as.data.frame(rlist[1]) %>% rownames_to_column("gene")
  r1 <- na.omit(r1[r1$padj < padj_threshold,])
  r2 <- as.data.frame(rlist[2]) %>% rownames_to_column("gene")
  r2 <- na.omit(r2[r2$padj < padj_threshold,])
  r3 <- as.data.frame(rlist[3]) %>% rownames_to_column("gene")
  r3 <- na.omit(r3[r3$padj < padj_threshold,])
  r <- rbind(r1,r2,r3)
  genes <- r$gene
  names(r)[names(r) == "gene"] <- "transcript"
  genes <- substr(genes, start = 1, stop = 9)
  r$gene <- genes

  subset_annotated <- biomartr::biomart(genes = r$gene,
                                        mart = "plants_mart",                 
                                        dataset = "athaliana_eg_gene",           
                                        attributes = attr,        
                                        filters =  "ensembl_gene_id" )
  names(subset_annotated)[names(subset_annotated) == 'ensembl_gene_id'] <- 'gene'
  
  res_annotated <- left_join(r, subset_annotated, by = "gene")
  res_annotated <- res_annotated[!(res_annotated$go_id == ""),]
  return(res_annotated)
}
