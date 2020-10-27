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
dds_2n <- DESeq(dds_2n)
res_2n <- results(dds_2n)
head(res_2n)
res_2n$genes <- rownames(res_2n)
res_2n_f <- as.data.frame(res_2n) %>% filter(padj < 0.05)
dim(res_2n_f)
hist(res_2n$padj, col="lightblue", main = "Adjusted p-value distribution")
diff_2n <- res_2n_f %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange), desc(padj))
write.table(diff_2n, file = "deg_2n_0.05.csv", sep=";", dec=".", row.names = TRUE)

#make volcano plot
resultsNames(dds_2n)
res_2n_shr <- lfcShrink(dds = dds_2n,
                        res = res_2n,
                        type = "apeglm",
                        coef = "treatment_n_vs_e")

EnhancedVolcano(toptable = res_2n_shr, 
                x = "log2FoldChange",           
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                lab = rownames(res_2n_shr),
                title = 'DEGs t=2 1.5nM solA treatment',
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

res_2n_f_annotated <- biomartr::biomart(genes = res_2n_f$genes,
                                        mart       = "plants_mart",                 
                                        dataset    = "athaliana_eg_gene",           
                                        attributes = attributes_to_retrieve,        
                                        filters =  "ensembl_gene_id" )
colnames(diff_2n) <- c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", "ensembl_gene_id")
diff_2n_annotated <- left_join(res_2n_f_annotated, diff_2n, by = "ensembl_gene_id")

write_delim(x = diff_2n_annotated,
            file = "DEGs_2n_annotated.csv",
            delim = ';')

#assign GO terms to DEGs
ora_2n_MF <- enrichGO(gene = res_2n_f_annotated$entrezgene_id, 
                       universe = all_arabidopsis_genes_annotated$entrezgene_id, 
                       OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                       keyType = "ENTREZID",
                       ont = "MF",              # either "BP", "CC" or "MF",
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       readable = TRUE, 
                       pool = FALSE)
ora_2n_MF@result[1:5,1:8]
write_delim(x = as.data.frame(ora_2n_MF@result), 
            file = "go_results_2n_MF.tsv", 
            delim = "\t")
dotplot(ora_2n_MF)

#KEGG
search_kegg_organism('ath', by = 'kegg_code')
search_kegg_organism('Arabidopsis thaliana', by = 'scientific_name')

ora_kegg <- enrichKEGG(gene = res_2n_f_annotated$entrezgene_id,
                       universe = all_arabidopsis_genes_annotated$entrezgene_id,
                       organism = "ath",
                       keyType = "ncbi-geneid",
                       minGSSize = 10,
                       maxGSSize = 500,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       use_internal_data = FALSE) # force to query latest KEGG db
dotplot(ora_kegg,
        color = "qvalue",
        showCategory = 10,
        size = "Count")

#KEGG modules
ora_kegg_mod <- enrichMKEGG(gene = res_2m_f_annotated$entrezgene_id,
                            universe = all_arabidopsis_genes_annotated$entrezgene_id,
                            organism = "ath",
                            keyType = "ncbi-geneid",
                            minGSSize = 10,           # minimal size of genes annotated by Ontology term for testing.
                            maxGSSize = 500,          # maximal size of genes annotated for testing
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05)
dotplot(ora_kegg_mod, 
        color = "p.adjust", 
        showCategory = 10, 
        size = "Count")

#iPath v3
diff_2m_annotated %>%
  filter(uniprotswissprot != "") %>%                                       # to remove genes with no matching Uniprot entries
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>%   # to create an ID that iPath can use
  dplyr::select(id_for_ipath) %>%                                          # we keep only the relevant ID for further copy-pasting 
  write.table(., 
              file = "diff_genes_swissprot.tsv", 
              row.names = FALSE, 
              quote = FALSE)
#throw resulting table into https://pathways.embl.de/