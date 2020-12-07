setwd("C:/Users/levla/github/Ath_RNAseq/Scheible_data")
dir()
library(limma)
treatment <- c(rep("FN", 3), rep("NoN", 3))
sample <- c("FN_A", "FN_B", "FN_C", "NoN_A", "NoN_B", "NoN_C")
filename <- c("FN_A.CEL", "FN_B.CEL", "FN_C.CEL", "NoN_A.CEL", "NoN_B.CEL", "NoN_C.CEL")
targets <- data.frame(sample, filename, treatment)
targets <- targets %>% column_to_rownames("sample")
library(affy)
eset <- justRMA(filenames = targets$filename)
colnames(eset) <- row.names(targets)
head(exprs(eset))
plotMDS(eset)
treatment <- factor(targets$treatment)
design <- model.matrix(~treatment)
fit <- lmFit(eset, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
results <- decideTests(fit)
summary(results)
plotMD(fit)
DEGs_scheible <- topTable(fit, n = 12000)
DEGs_scheible <- DEGs_scheible %>% dplyr::filter(adj.P.Val < 0.05)
DEGs_scheible <- DEGs_scheible %>% rownames_to_column("Affy_identifier")
all_arabidopsis_genes_annotated <- get_Athaliana_annotations(attr = c("affy_ath1_121501", "entrezgene_id"))
names(all_arabidopsis_genes_annotated)[2] <- "Affy_identifier"
combined <- left_join(DEGs_scheible, all_arabidopsis_genes_annotated, by = "Affy_identifier")
combined <- combined[,-c(3,4,5,7,9)]
names(combined) <- c("Affy_identifier", "logFC_scheible2004", "Padj_scheible2004", "genes")
degs <- left_join(alldegs, combined, by = "genes")

#allDEGs is the 'clean' DEG list of my experiment
#degs is the fusion between scheible and my DEG list
#combined are the annotated DEGs from scheible only
scheibl_in_own <- combined[combined$genes %in% degs$genes,]
names(all_arabidopsis_genes_annotated)[1] <- "genes"
scheibl_in_own <- left_join(scheibl_in_own, all_arabidopsis_genes_annotated, by = "genes")

go <- enrichGO(gene = scheibl_in_own$entrezgene_id, 
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
write.csv(degs, "degs.csv")
l <- unique(degs$genes) == unique(allDEGs$genes)
