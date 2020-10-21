library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(DESeq2)
library(PCAtools)
#make normalised dataset + plot PCA
dds_norm <- estimateSizeFactors(dds)
dds_norm_vst <- vst(dds_norm)
counts_norm_vst <- assay(dds_norm_vst, normalized = TRUE)

#do prcomp with only top 500 variable genes
rv <- rowVars(assay(dds_norm_vst))
ntop <- 500
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
counts_norm_vst_t_500 <- t( assay(dds_norm_vst)[select, ] )
dds_norm_vst_prc500 <-prcomp(counts_norm_vst_t_500)
dds_norm_vst_prc500_df <- as.data.frame(dds_norm_vst_prc500$x)
dds_norm_vst_prc500_df$PC1
pcaData <- plotPCA(dds_norm_vst, intgroup = c("treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(dds_norm_vst_prc500_df, aes(PC1, PC2, color = time, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

#make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes)
pcaData <- plotPCA(dds_norm_vst, intgroup = c("treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = time)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

#do PCA analysis with only top 500 variable genes with PCAtools
xp_design_pca <- xp_design
row.names(xp_design_pca) <- xp_design[,1]
xp_design_pca <- xp_design_pca[,-1]
remvar <- ((nrow(counts) - ntop)/nrow(counts))
dds_vst_pca <- pca(counts_norm_vst, metadata = xp_design_pca, removeVar = remvar)

#make scree plot (500 top genes)
percentVar_all <- dds_norm_vst_prc500$sdev^2 / sum(dds_norm_vst_prc500$sdev^2)
scree_all <- data.frame(percentVar_all)
scree_all[,2] <- c(1:46)
colnames(scree_all) <- c("variance","component_nr")
ggplot(scree_all, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for all timepoints")


#make bi plot (500 top) with PCAtools
biplot(dds_vst_pca,
       colby = "time",
       shape = "treatment")

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca$variance)

#make pairs plot (500 top)
pairsplot(dds_vst_pca,
          components = getComponents(dds_vst_pca, c(1:6)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))


#-------------------------------------------------------------------------------------------------------------
#do same type of analsysi as above but with 50% of the highest variance genes

#do prcomp with only top 50% variable genes
ntop_.5 <- nrow(counts)/2
ntop_.5 <- as.integer(ntop_.5)
select_.5 <- order(rv, decreasing = TRUE)[seq_len(min(ntop_.5, length(rv)))]
counts_norm_vst_t_.5 <- t(assay(dds_norm_vst)[select_.5, ] )
dds_norm_vst_prc_.5 <-prcomp(counts_norm_vst_t_.5)
dds_norm_vst_prc_.5_df <- as.data.frame(dds_norm_vst_prc_.5$x)
percentVar_all_.5 <- dds_norm_vst_prc_.5$sdev^2 / sum(dds_norm_vst_prc_.5$sdev^2)
ggplot(dds_norm_vst_prc_.5_df, aes(PC1, PC2, color = time)) +
  geom_point() +
  xlab(paste0("PC1: ",percentVar_all_.5[1]*100,"% variance")) +
  ylab(paste0("PC2: ",percentVar_all_.5[2]*100,"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with only top 50% variable genes with PCAtools
dds_vst_pca_.5 <- pca(counts_norm_vst, metadata = xp_design_pca, removeVar = .5)

#make scree plot (50% top genes)
scree_all_.5 <- data.frame(percentVar_all_.5)
scree_all_.5[,2] <- c(1:46)
colnames(scree_all_.5) <- c("variance","component_nr")
ggplot(scree_all_.5, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for all timepoints")


#make bi plot (50% top) with PCAtools
biplot(dds_vst_pca_.5,
       colby = "time",
       shape = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_.5$variance)

#make pairs plot (50% top)
pairsplot(dds_vst_pca_.5,
          components = getComponents(dds_vst_pca_.5, c(1:10)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#---------------------------------------------------------------------------------------------------------------
#do same analyses but only removing the lowest 10% of variation (which is default of PCAtools package)

#do prcomp with lower 10% of variable genes removed
ntop_.1 <- nrow(counts)*.9
ntop_.1 <- as.integer(ntop_.1)
select_.1 <- order(rv, decreasing = TRUE)[seq_len(min(ntop_.1, length(rv)))]
counts_norm_vst_t_.1 <- t( assay(dds_norm_vst)[select_.1, ] )
dds_norm_vst_prc_.1 <-prcomp(counts_norm_vst_t_.1)
dds_norm_vst_prc_.1_df <- as.data.frame(dds_norm_vst_prc_.1$x)
percentVar_all_.1 <- dds_norm_vst_prc_.1$sdev^2 / sum(dds_norm_vst_prc_.1$sdev^2)
ggplot(dds_norm_vst_prc_.1_df, aes(PC1, PC2, color = time, shape = treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentVar_all_.1[1]*100,"% variance")) +
  ylab(paste0("PC2: ",percentVar_all_.1[2]*100,"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis, removing lower 10% of variance
dds_vst_pca_.1 <- pca(counts_norm_vst, metadata = xp_design_pca, removeVar = .1)

#make scree plot (90% top genes)
scree_all_.1 <- data.frame(percentVar_all_.1)
scree_all_.1[,2] <- c(1:46)
colnames(scree_all_.1) <- c("variance","component_nr")
ggplot(scree_all_.1, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for all timepoints")


#make bi plot (90% top) with PCAtools
biplot(dds_vst_pca_.1,
       colby = "time",
       shape = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_.1$variance)

#make pairs plot (90% top)
pairsplot(dds_vst_pca_.1,
          components = getComponents(dds_vst_pca_.1, c(1:10)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#-------------------------------------------------------------------------------------------------------
#above scripts did unscaled PCA, try also scaled PCA (with PCAtools)

#do prcomp with lower 10% of variable genes removed
select_.1 <- order(rv, decreasing = TRUE)[seq_len(min(ntop_.1, length(rv)))]
counts_norm_vst_t_.1 <- t( assay(dds_norm_vst)[select_.1, ] )
dds_norm_vst_prc_.1 <-prcomp(counts_norm_vst_t_.1, scale = TRUE)
dds_norm_vst_prc_.1_df <- as.data.frame(dds_norm_vst_prc_.1$x)
percentVar_all_.1 <- dds_norm_vst_prc_.1$sdev^2 / sum(dds_norm_vst_prc_.1$sdev^2)
ggplot(dds_norm_vst_prc_.1_df, aes(PC1, PC2, color = time, shape = treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentVar_all_.1[1]*100,"% variance")) +
  ylab(paste0("PC2: ",percentVar_all_.1[2]*100,"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with PCAtools, removing lower 10% of variance
ntop_.1 <- nrow(counts)*.9
ntop_.1 <- as.integer(ntop_.1)
dds_vst_pca_.1 <- pca(counts_norm_vst, metadata = xp_design_pca, removeVar = .1, scale = TRUE, center = TRUE)

#make scree plot (90% top genes)
scree_all_.1 <- data.frame(percentVar_all_.1)
scree_all_.1[,2] <- c(1:46)
colnames(scree_all_.1) <- c("variance","component_nr")
ggplot(scree_all_.1, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for all timepoints")


#make bi plot (90% top) with PCAtools
biplot(dds_vst_pca_.1,
       colby = "time",
       shape = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_.1$variance)

#make pairs plot (90% top)
pairsplot(dds_vst_pca_.1,
          components = getComponents(dds_vst_pca_.1, c(1:9)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#-------------------------------------------------------------------------------------------------------
#take only 500 top variance genes but this time scaled

#do prcomp with only top 500 variable genes
ntop <- 500
rv <- rowVars(assay(dds_norm_vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
counts_norm_vst_t_500 <- t(assay(dds_norm_vst)[select, ] )
dds_norm_vst_prc500 <- prcomp(counts_norm_vst_t_500, scale = TRUE)
dds_norm_vst_prc500_df <- as.data.frame(dds_norm_vst_prc500$x)
percentVar_all <- dds_norm_vst_prc500$sdev^2 / sum(dds_norm_vst_prc500$sdev^2)
ggplot(dds_norm_vst_prc500_df, aes(PC1, PC2, color = time, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_all[1]*100,"% variance")) +
  ylab(paste0("PC2: ",percentVar_all[2]*100,"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

#do PCA analysis with only top 500 variable genes with PCAtools
xp_design_pca <- xp_design
row.names(xp_design_pca) <- xp_design[,1]
xp_design_pca <- xp_design_pca[,-1]
remvar <- ((nrow(counts) - ntop)/nrow(counts))
dds_vst_pca <- pca(counts_norm_vst, metadata = xp_design_pca, removeVar = remvar, scale = TRUE)

#make scree plot (500 top genes)
scree_all <- data.frame(percentVar_all)
scree_all[,2] <- c(1:46)
colnames(scree_all) <- c("variance","component_nr")
ggplot(scree_all, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for all timepoints")


#make bi plot (500 top) with PCAtools
biplot(dds_vst_pca,
       colby = "time",
       shape = "treatment")

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca$variance)

#make pairs plot (500 top)
pairsplot(dds_vst_pca,
          components = getComponents(dds_vst_pca, c(1:6)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))
