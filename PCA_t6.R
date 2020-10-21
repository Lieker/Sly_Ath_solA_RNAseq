library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(DESeq2)
library(PCAtools)
#make normalised dataset + plot PCA
dds_norm_6 <- estimateSizeFactors(dds_6)
dds_norm_vst_6 <- vst(dds_norm_6)
counts_norm_vst_6 <- assay(dds_norm_vst_6, normalized = TRUE)

#do prcomp with only top 500 variable genes
rv6 <- rowVars(assay(dds_norm_vst_6))
ntop <- 500
select6 <- order(rv6, decreasing = TRUE)[seq_len(min(ntop, length(rv6)))]
counts_norm_vst_t6_500 <- t( assay(dds_norm_vst_6)[select6, ] )
dds_norm_vst_6_prc500 <-prcomp(counts_norm_vst_t6_500)
dds_norm_vst_6_prc500_df <- as.data.frame(dds_norm_vst_6_prc500$x)
percentVar_6 <- dds_norm_vst_6_prc500$sdev^2 / sum(dds_norm_vst_6_prc500$sdev^2) *100
dds_norm_vst_6_prc500_df$treatment <- rep(c("e","m","n"), c(5,5,5))
ggplot(dds_norm_vst_6_prc500_df, aes(PC1, PC2, color=treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_6[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_6[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot for t=6")

#make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes)
pcaData <- plotPCA(dds_norm_vst_6, intgroup = c("treatment"), returnData = TRUE)
percentVar_6_plotpca <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_6_plotpca[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_6_plotpca[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

#do PCA analysis with only top 500 variable genes with PCAtools
xp_design_6_pca <- xp_design_6
row.names(xp_design_6_pca) <- xp_design_6[,1]
xp_design_6_pca <- xp_design_6_pca[,-1]
remvar <- ((nrow(counts_6) - ntop)/nrow(counts_6))
dds_vst_pca_6 <- pca(counts_norm_vst_6, metadata = xp_design_6_pca, removeVar = remvar)

#make scree plot (500 top genes)
scree_6 <- data.frame(percentVar_6)
scree_6[,2] <- c(1:15)
colnames(scree_6) <- c("variance","component_nr")
ggplot(scree_6, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=6")


#make bi plot (500 top) with PCAtools
biplot(dds_vst_pca_6,
       colby = "treatment")

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_6$variance)

#make pairs plot (500 top)
pairsplot(dds_vst_pca_6,
          components = getComponents(dds_vst_pca_6, c(1:5)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))


#-------------------------------------------------------------------------------------------------------------
#do same type of analsysi as above but with 50% of the highest variance genes

#do prcomp with only top 50% variable genes
ntop_.5 <- nrow(counts_6)/2
ntop_.5 <- as.integer(ntop_.5)
select_6_.5 <- order(rv6, decreasing = TRUE)[seq_len(min(ntop_.5, length(rv6)))]
counts_norm_vst_t6_.5 <- t(assay(dds_norm_vst_6)[select_6_.5, ] )
dds_norm_vst_6_prc_.5 <-prcomp(counts_norm_vst_t6_.5)
dds_norm_vst_6_prc_.5_df <- as.data.frame(dds_norm_vst_6_prc_.5$x)
percentVar_6_.5 <- dds_norm_vst_6_prc_.5$sdev^2 / sum(dds_norm_vst_6_prc_.5$sdev^2) *100
dds_norm_vst_6_prc_.5_df$treatment <- rep(c("e","m","n"), c(5,5,5))
ggplot(dds_norm_vst_6_prc_.5_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_6_.5[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_6_.5[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with only top 50% variable genes with PCAtools
dds_vst_pca_6_.5 <- pca(counts_norm_vst_6, metadata = xp_design_6_pca, removeVar = .5)

#make scree plot (50% top genes)
scree_6_.5 <- data.frame(percentVar_6_.5)
scree_6_.5[,2] <- c(1:15)
colnames(scree_6_.5) <- c("variance","component_nr")
ggplot(scree_6_.5, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=6")


#make bi plot (50% top) with PCAtools
biplot(dds_vst_pca_6_.5,
       colby = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_6_.5$variance)

#make pairs plot (50% top)
pairsplot(dds_vst_pca_6_.5,
          components = getComponents(dds_vst_pca_6_.5, c(1:4)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#-------------------------------------------------------------------------------------------------------------------
#do scaled analysis with only lowest 10% removed

#do prcomp with lower 10% of variable genes removed
ntop_.1 <- nrow(counts)*.9
ntop_.1 <- as.integer(ntop_.1)
select_6_.1 <- order(rv6, decreasing = TRUE)[seq_len(min(ntop_.1, length(rv6)))]
counts_norm_vst_t6_.1 <- t( assay(dds_norm_vst_6)[select_6_.1, ] )
dds_norm_vst_6_prc_.1 <-prcomp(counts_norm_vst_t6_.1, scale = TRUE)
dds_norm_vst_6_prc_.1_df <- as.data.frame(dds_norm_vst_6_prc_.1$x)
percentVar_6_.1 <- dds_norm_vst_6_prc_.1$sdev^2 / sum(dds_norm_vst_6_prc_.1$sdev^2) *100
dds_norm_vst_6_prc_.1_df$treatment <- rep(c("e","m","n"), c(5,5,5))
ggplot(dds_norm_vst_6_prc_.1_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_6_.1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_6_.1[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with PCAtools, removing lower 10% of variance
dds_vst_pca_6_.1 <- pca(counts_norm_vst_6, metadata = xp_design_6_pca, removeVar = .1, scale = TRUE, center = TRUE)

#make scree plot (90% top genes)
scree_6_.1 <- data.frame(percentVar_6_.1)
scree_6_.1[,2] <- c(1:15)
colnames(scree_6_.1) <- c("variance","component_nr")
ggplot(scree_6_.1, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=6")


#make bi plot (90% top) with PCAtools
biplot(dds_vst_pca_6_.1,
       colby = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_6_.1$variance)

#make pairs plot (90% top)
pairsplot(dds_vst_pca_6_.1,
          components = getComponents(dds_vst_pca_6_.1, c(1:3)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#--------------------------------------------------------------------------------------------------------

#permanova for t=6 subdataset
library(vegan)
counts_norm_vst_6_t <- as.data.frame(t(counts_norm_vst_6))
adonis(counts_norm_vst_6_t~treatment, data=xp_design_6, permutations = 9999, method = "euclidean")
