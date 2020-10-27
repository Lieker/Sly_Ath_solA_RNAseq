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
dds_24m <- DESeq(dds_24m)
res_24m <- results(dds_24m)
head(res_24m)
res_24m$genes <- rownames(res_24m)
res_24m_f <- as.data.frame(res_24m) %>% filter(padj < 0.05)
dim(res_24m_f)
hist(res_24m$padj, col="lightblue", main = "Adjusted p-value distribution")
diff_24m <- res_24m_f %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange), desc(padj))
write.table(diff_24m, file = "deg_24m_0.05.csv", sep=";", dec=".", row.names = TRUE)

#make volcano plot
resultsNames(dds_24m)
res_24m_shr <- lfcShrink(dds = dds_24m,
                        res = res_24m,
                        type = "apeglm",
                        coef = "treatment_m_vs_e")

EnhancedVolcano(toptable = res_24m_shr, 
                x = "log2FoldChange",           
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                lab = rownames(res_24m_shr),
                title = 'DEGs t=24 1.5uM solA treatment',
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

res_24m_f_annotated <- biomartr::biomart(genes = res_24m_f$genes,
                                        mart       = "plants_mart",                 
                                        dataset    = "athaliana_eg_gene",           
                                        attributes = attributes_to_retrieve,        
                                        filters =  "ensembl_gene_id" )
colnames(diff_24m) <- c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", "ensembl_gene_id")
diff_24m_annotated <- left_join(res_24m_f_annotated, diff_24m, by = "ensembl_gene_id")

write_delim(x = diff_24m_annotated,
            file = "DEGs_24m_annotated.csv",
            delim = ';')

#assign GO terms to DEGs
ora_24m_MF <- enrichGO(gene = res_24m_f_annotated$entrezgene_id, 
                       universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                       OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                       keyType = "ENTREZID",
                       ont = "MF",              # either "BP", "CC" or "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE, 
                       pool = FALSE)
ora_24m_MF@result[1:5,1:8]
write_delim(x = as.data.frame(ora_24m_MF@result), 
            file = "go_results_24m_MF.tsv", 
            delim = "\t")
dotplot(ora_24m_MF)
