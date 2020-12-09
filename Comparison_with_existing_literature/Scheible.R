library(dplyr)
library(tidyr)
library(tidyverse)
library(affy)
library(ecolicdf)
library(biomartr)
library("biomaRt")
library(org.At.tair.db)
library(limma)
library(clusterProfiler)
library(org.At.tair.db)
scheible <- read.csv("Scheibletable2.csv", header = TRUE, sep = ",")
names(scheible)[1] <- "Affy_identifier"
s <- scheible[,1:13]
s <- s %>% column_to_rownames(var = "Affy_identifier")
biomartr::organismBM(organism = "Arabidopsis thaliana")
attr = biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")
all_arabidopsis_genes <- read.csv("raw_counts.csv", header = TRUE, stringsAsFactors = FALSE)
all_arabidopsis_genes <- all_arabidopsis_genes$Geneid
all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                     mart       = "plants_mart",                 
                                                     dataset    = "athaliana_eg_gene",           
                                                     attributes = c("affy_ath1_121501",
                                                                    "entrezgene_id"),        
                                                     filters =  "ensembl_gene_id" ) 
names(s) <- c("T0_noadd",
              "30min_noadd",
              "3h_noadd",
              "30min_KCl",
              "3h_KCl",
              "30min_KNO3_A",
              "30min_KNO3_B",
              "3h_KNO3_A",
              "3h_KNO3_B",
              "fullnut_A",
              "fullnut_B",
              "fullnut_C")
time <- c(0,30,180,30,180,30,30,180,180,0,0,0)
treatment <- c("no","no","no","KCl","KCl","KNO3","KNO3","KNO3","KNO3","no","no","no")
nutrition <- c("noN","noN","noN","noN","noN","noN","noN","noN","noN","yesN","yesN","yesN")
xp_design <- data.frame(names(s),time,treatment,nutrition)
write.csv(xp_design, "xp_design_Scheible.csv", row.names = FALSE)
targets <- readTargets("xp_design_Scheible.csv", sep = ",")

eset <- s
head(eset)
plotMDS(eset)
targets_1 <- targets %>% dplyr::filter(treatment == "no")
eset_1 <- eset[,-c(4:9)]
nutrition_1 <- factor(targets_1$nutrition)
design <- model.matrix(~nutrition_1)
fit <- lmFit(eset_1, design)
fit <-  eBayes(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit)
summary(results)
degs_s <- topTable(fit, n = 7000)
degs_s <- degs_s %>% dplyr::filter(degs_s$adj.P.Val < 0.05 )
degs_s <- degs_s %>% rownames_to_column()
names(degs_s)[1] <- "Affy_identifier"
names(all_arabidopsis_genes_annotated) <- c("ensembl_geneid", "Affy_identifier","entrezgene_id")
degs_s <- inner_join(degs_s, all_arabidopsis_genes_annotated, by = "Affy_identifier")
write.csv(degs_s, "Scheible_DEGs.csv", row.names = FALSE)
degs <- read.csv("allDEGs.csv", header = TRUE, sep = ";")
combined_with_scheible <- degs[degs$genes %in% degs_s$ensembl_geneid,]
lfc_scheible <- degs_s[degs_s$ensembl_geneid %in% degs$genes,]
names(lfc_scheible)
names(lfc_scheible)[2] <- "logFC_scheible"
lfc_scheible <- lfc_scheible[,c(1,2,6,8,9)]
names(lfc_scheible)[4] <- "genes"
combined_with_scheible <- inner_join(combined_with_scheible, lfc_scheible, by = "genes")
names(combined_with_scheible)[17] <- "adj_Pval_scheible"
write.csv(combined_with_scheible, "Scheible_combined_with_own.csv", row.names = FALSE)
go <- enrichGO(gene = combined_with_scheible$entrezgene_id.x, 
               universe = all_arabidopsis_genes_annotated$entrezgene_id, 
               OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
               keyType = "ENTREZID",
               ont = "BP",              # either "BP", "CC" or "MF",
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05,
               readable = TRUE, 
               pool = FALSE)
goo <- go@result
write.csv(goo, "go_combi_scheible.csv")
