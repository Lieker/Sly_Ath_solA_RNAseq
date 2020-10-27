BiocManager::install("EnhancedVolcano")
BiocManager::install("apeglm")
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("biomaRt")
install.packages("biomartr", dependencies = TRUE)
install.packages("clusterProfiler")
BiocManager::install("org.At.tair.db")
BiocManager::install("clusterProfiler")

library(apeglm)
library(EnhancedVolcano)
library(biomartr)
library(tidyverse)
library("biomaRt")
library(clusterProfiler)
library(org.At.tair.db)

#do DEG analysis
dds_24n <- DESeq(dds_24n)
res_24n <- results(dds_24n)
head(res_24n)
res_24n$genes <- rownames(res_24n)
res_24n_f <- as.data.frame(res_24n) %>% filter(padj < 0.05)
dim(res_24n_f)
hist(res_24n$padj, col="lightblue", main = "Adjusted p-value distribution")
diff_24n <- res_24n_f %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange), desc(padj))
write.table(diff_24n, file = "deg_24n_0.05.csv", sep=";", dec=".", row.names = TRUE)

#make volcano plot
resultsNames(dds_24n)
res_24n_shr <- lfcShrink(dds = dds_24n,
                        res = res_24n,
                        type = "apeglm",
                        coef = "treatment_n_vs_e")

EnhancedVolcano(toptable = res_24n_shr, 
                x = "log2FoldChange",           
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                lab = rownames(res_24n_shr),
                title = 'DEGs t=24 1.5nM solA treatment',
                subtitle = 'LFC>1, padj <0.05')

#add annotations to DEGs
biomartr::organismBM(organism = "Arabidopsis thaliana")
arabido_attributes = 
  biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")
arabido_attributes

attributes_to_retrieve <- c("tair_symbol", 
                            "uniprotswissprot",
                            "entrezgene_id",
                            "description",
                            "external_gene_name",
                            "external_gene_source")
all_arabidopsis_genes <- rownames(counts)
all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                     mart       = "plants_mart",                 
                                                     dataset    = "athaliana_eg_gene",           
                                                     attributes = attributes_to_retrieve,        
                                                     filters =  "ensembl_gene_id" )  
all_arabidopsis_genes_annotated$entrezgene_id = as.character(
  all_arabidopsis_genes_annotated$entrezgene_id) 

res_24n_f_annotated <- biomartr::biomart(genes = res_24n_f$genes,
                                        mart       = "plants_mart",                 
                                        dataset    = "athaliana_eg_gene",           
                                        attributes = attributes_to_retrieve,        
                                        filters =  "ensembl_gene_id" )
colnames(diff_24n) <- c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", "ensembl_gene_id")
diff_24n_annotated <- left_join(res_24n_f_annotated, diff_24n, by = "ensembl_gene_id")

write_delim(x = diff_24n_annotated,
            file = "DEGs_24n_annotated.csv",
            delim = ';')

#assign GO terms to DEGs
ora_24n_MF <- enrichGO(gene = res_24n_f_annotated$entrezgene_id, 
                       universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                       OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                       keyType = "ENTREZID",
                       ont = "MF",              # either "BP", "CC" or "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE, 
                       pool = FALSE)
ora_24n_MF@result[1:5,1:8]
write_delim(x = as.data.frame(ora_24n_MF@result), 
            file = "go_results_24n_MF.tsv", 
            delim = "\t")
dotplot(ora_24n_MF)
