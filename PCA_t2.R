library(ggplot2)
library(tidyr)
library(dplyr)
library(matrixStats)
library(DESeq2)
library(PCAtools)
#make normalised dataset + plot PCA
dds_norm_2 <- estimateSizeFactors(dds_2)
dds_norm_vst_2 <- vst(dds_norm_2)
counts_norm_vst_2 <- assay(dds_norm_vst_2, normalized = TRUE)

#do prcomp with only top 500 variable genes
rv2 <- rowVars(assay(dds_norm_vst_2))
ntop <- 500
select2 <- order(rv2, decreasing = TRUE)[seq_len(min(ntop, length(rv2)))]
counts_norm_vst_t2_500 <- t( assay(dds_norm_vst_2)[select2, ] )
dds_norm_vst_2_prc500 <-prcomp(counts_norm_vst_t2_500)
dds_norm_vst_2_prc500_df <- as.data.frame(dds_norm_vst_2_prc500$x)
percentVar_2 <- dds_norm_vst_2_prc500$sdev^2 / sum(dds_norm_vst_2_prc500$sdev^2) *100
dds_norm_vst_2_prc500_df$treatment <- rep(c("e","m","n"), c(6,5,5))
ggplot(dds_norm_vst_2_prc500_df, aes(PC1, PC2, color=treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_2[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot for t=2")

#make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes)
pcaData <- plotPCA(dds_norm_vst_2, intgroup = c("treatment"), returnData = TRUE)
percentVar_2_plotpca <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_2_plotpca[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_2_plotpca[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

#do PCA analysis with only top 500 variable genes with PCAtools
xp_design_2_pca <- xp_design_2
row.names(xp_design_2_pca) <- xp_design_2[,1]
xp_design_2_pca <- xp_design_2_pca[,-1]
remvar <- ((nrow(counts_2) - ntop)/nrow(counts_2))
dds_vst_pca_2 <- pca(counts_norm_vst_2, metadata = xp_design_2_pca, removeVar = remvar)

#make scree plot (500 top genes)
scree_2 <- data.frame(percentVar_2)
scree_2[,2] <- c(1:16)
colnames(scree_2) <- c("variance","component_nr")
ggplot(scree_2, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=2")


#make bi plot (500 top) with PCAtools
biplot(dds_vst_pca_2,
       colby = "treatment")

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_2$variance)

#make pairs plot (500 top)
pairsplot(dds_vst_pca_2,
          components = getComponents(dds_vst_pca_2, c(1:6)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))


#-------------------------------------------------------------------------------------------------------------
#do same type of analsysi as above but with 50% of the highest variance genes

#do prcomp with only top 50% variable genes
ntop_.5 <- nrow(counts_2)/2
ntop_.5 <- as.integer(ntop_.5)
select_2_.5 <- order(rv2, decreasing = TRUE)[seq_len(min(ntop_.5, length(rv2)))]
counts_norm_vst_t2_.5 <- t(assay(dds_norm_vst_2)[select_2_.5, ] )
dds_norm_vst_2_prc_.5 <-prcomp(counts_norm_vst_t2_.5)
dds_norm_vst_2_prc_.5_df <- as.data.frame(dds_norm_vst_2_prc_.5$x)
percentVar_2_.5 <- dds_norm_vst_2_prc_.5$sdev^2 / sum(dds_norm_vst_2_prc_.5$sdev^2) *100
dds_norm_vst_2_prc_.5_df$treatment <- rep(c("e","m","n"), c(6,5,5))
ggplot(dds_norm_vst_2_prc_.5_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar_2_.5[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_2_.5[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with only top 50% variable genes with PCAtools
dds_vst_pca_2_.5 <- pca(counts_norm_vst_2, metadata = xp_design_2_pca, removeVar = .5)

#make scree plot (50% top genes)
scree_2_.5 <- data.frame(percentVar_2_.5)
scree_2_.5[,2] <- c(1:16)
colnames(scree_2_.5) <- c("variance","component_nr")
ggplot(scree_2_.5, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=2")


#make bi plot (50% top) with PCAtools
biplot(dds_vst_pca_.5,
       colby = "time",
       shape = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_2_.5$variance)

#make pairs plot (50% top)
pairsplot(dds_vst_pca_2_.5,
          components = getComponents(dds_vst_pca_2_.5, c(1:5)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#-------------------------------------------------------------------------------------------------------------------
#do scaled analysis with only lowest 10% removed

#do prcomp with lower 10% of variable genes removed
ntop_.1 <- nrow(counts)*.9
ntop_.1 <- as.integer(ntop_.1)
select_2_.1 <- order(rv2, decreasing = TRUE)[seq_len(min(ntop_.1, length(rv2)))]
counts_norm_vst_t2_.1 <- t( assay(dds_norm_vst_2)[select_2_.1, ] )
dds_norm_vst_2_prc_.1 <-prcomp(counts_norm_vst_t2_.1, scale = TRUE)
dds_norm_vst_2_prc_.1_df <- as.data.frame(dds_norm_vst_2_prc_.1$x)
percentVar_2_.1 <- dds_norm_vst_2_prc_.1$sdev^2 / sum(dds_norm_vst_2_prc_.1$sdev^2) *100
dds_norm_vst_2_prc_.1_df$treatment <- rep(c("e","m","n"), c(6,5,5))
ggplot(dds_norm_vst_2_prc_.1_df, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentVar_2_.1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_2_.1[2],"% variance")) + 
  coord_fixed() +
  ggtitle("PCA plot")

##make bi plot acc to RNAseq lesson (same as prcomp 500 most variable genes) is not possible, because default is 500 genes

#do PCA analysis with PCAtools, removing lower 10% of variance
dds_vst_pca_2_.1 <- pca(counts_norm_vst_2, metadata = xp_design_2_pca, removeVar = .1, scale = TRUE, center = TRUE)

#make scree plot (90% top genes)
scree_2_.1 <- data.frame(percentVar_2_.1)
scree_2_.1[,2] <- c(1:16)
colnames(scree_2_.1) <- c("variance","component_nr")
ggplot(scree_2_.1, mapping=aes(x=component_nr, y=variance)) + 
  geom_bar(stat="identity") +
  ggtitle("screeplot for t=2")


#make bi plot (90% top) with PCAtools
biplot(dds_vst_pca_2_.1,
       colby = "treatment",
       pointSize = 5,
       legendPosition = 'right')

#find PC's to keep
elbow <- findElbowPoint(dds_vst_pca_2_.1$variance)

#make pairs plot (90% top)
pairsplot(dds_vst_pca_2_.1,
          components = getComponents(dds_vst_pca_2_.1, c(1:5)),
          colby = "treatment",
          colkey = c("e" = "gray",
                     "m" = "darkgreen",
                     "n" = "lightgreen"))

#----------------------------------------------------------------------------------------------------------------------
#since PC4 seems to show the division of treatment 'e' versus 'm' and 'n', we want to know which genes are
#responsible for this division --> loadingsplot

plotloadings(dds_vst_pca_2_.1,
             rangeRetain = 0.001,
             labSize = 3.0,
             title = 'Loadings plot t=2',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 0.1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

#try permanova on t=2 data
library(vegan)
counts_norm_vst_2_t <- as.data.frame(t(counts_norm_vst_2))
adonis(counts_norm_vst_2_t~treatment, data=xp_design_2, permutations = 9999, method = "euclidean")

#plsda on t=2 data
library(mixOmics)
Y = xp_design_2$treatment
summary(Y)
X = counts_norm_vst_2_t
head(X)
plsda <- plsda(X, Y, ncomp = 8)
background = background.predict(plsda, comp.predicted=2, dist = "max.dist", xlim = c(-150,200), ylim = c(-100,150)) 
plotIndiv(plsda,
          comp = 1:2,
          group = Y,
          ind.names = FALSE,
          legend = TRUE,
          ellipse = TRUE,
          title = "PLS-DA on t=2 data",
          background = background)
perf.plsda <- perf(plsda, validation = "Mfold", folds = 5, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 10) #takes a few minutes
plot(perf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
perf.plsda$choice.ncomp #tells you how many components you should keep in plsda(x,y,ncomp=)
#repeat plsda with ideal number of components
plsda <- plsda(X, Y, ncomp = 2)
plotIndiv(plsda,
          comp = 1:2,
          group = Y,
          ind.names = FALSE,
          legend = TRUE,
          ellipse = TRUE,
          title = "PLS-DA on t=2 data",
          background = background)
#now try sPLS-DA (only keeping a few genes that explain the treatments best)
list.keepX <- c(1:10, seq(20,300,10))
tune.splsda <- tune.splsda(X,Y,ncomp = 10, validation = 'Mfold', folds = 5,
                           progressBar = TRUE, dist = 'max.dist', measure = 'BER',
                           test.keepX = list.keepX, nrepeat = 10, cpus = 2)
                            #this will take about 10 min
error <- tune.splsda$error.rate
ncomp <- tune.splsda$choice.ncomp$ncomp
select.keepX <- tune.splsda$choice.keepX[1:ncomp]
select.keepX
plot(tune.splsda)
splsda <- splsda(X,Y,ncomp = 2, keepX = select.keepX)
plotIndiv(splsda, comp = c(1,2),
          group = Y, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on t=2')
auc.splsda = auroc(splsda, roc.comp = ncomp) #means that 2 components lead to perfect separation (100%)
perf <- perf(splsda, validation = 'Mfold', folds = 5,
             dist = 'max.dist', nrepeat = 10,
             progressBar = TRUE)
perf$error.rate
plot(perf)
#additional analysis plots
par(mfrow=c(1,2))
plot(perf$features$stable[[1]], type = 'h', ylab = 'stability',
     xlab = 'features', main = 'comp 1', las = 2)
plot(perf$features$stable[[2]], type = 'h', ylab = 'stability',
     xlab = 'features', main = 'comp 2', las = 2)
par(mfrow=c(1,1))

ind.match = match(selectVar(splsda, comp = 1)$name,
                  names(perf$features$stable[[1]]))
Freq = as.numeric(perf$features$stable[[1]][ind.match])
data.frame(selectVar(splsda, comp = 1)$value, Freq)
plotLoadings(splsda, comp = 1, title = 'Loadings on comp 1',
             contrib = 'max', method = 'mean')
cim(splsda)
cim(splsda, comp = 1, title = 'comp 1')
cim(splsda, comp = 2, title = 'comp 2')

plotArrow(splsda, legend = TRUE)

plotVar(splsda, overlap = F, cex = 2)
