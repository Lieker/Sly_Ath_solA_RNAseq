library(dplyr)
library(stringr)
library(DESeq2)

#assemble counts table and xp_design
setwd("Z:/Papers/Arabidopsis_RNAseq/R")
counts_E02_1_raw <- read.csv("counts_E02_1.csv", header = F, stringsAsFactors = F)
counts_E02_1 <- counts_E02_1_raw[-1,]
genes <- counts_E02_1[2:28776,1]
headers_E02_1 <- counts_E02_1[1,]
counts_E02_1 <- counts_E02_1[-1,]
counts_E02_1 <- counts_E02_1[,-1]
row.names(counts_E02_1) <- genes
counts_E02_1 <- counts_E02_1[,-c(1,2,3,4,5)]
names(counts_E02_1) <- c("E02_1_1", "E02_1_2")
#from here, you need to first run the code of the DEG_analysis script until (not including) #make dds dataset,
#in order to create a full dataset
counts_E <- counts[,-1]
counts_E$E02_1_1 <- counts_E02_1$E02_1_1
counts_E$E02_1_2 <- counts_E02_1$E02_1_2
batches_2E <- read.csv("batches_2E.csv", header = T, stringsAsFactors = F)
treatment_2E <- rep(c("e","m","n"), c(17, 15, 15))
time_2E <- rep(c("2","6","24","2","6","24","2","6","24"), c(7,5,5,5,5,5,5,5,5))
headers_E <- names(counts_E)
headers_E <- headers_E[order(headers_E)]
xp_design_2E <- data.frame(headers_E, treatment_2E, time_2E, batches_2E$batch, batches_2E$batch2, batches_2E$batch3)
colnames(xp_design_2E) <- c("sample","treatment","time","batch1","batch2","batch3")
all(colnames(counts_E) %in% xp_design_2E$sample)
counts_E <- counts_E[,order(colnames(counts_E))]
all(colnames(counts_E) == xp_design_2E$sample)
xp_design_2E$batch1 <- as.factor(xp_design_2E$batch1)
xp_design_2E$treatment <- as.factor(xp_design_2E$treatment)
xp_design_2E$batch3 <- as.factor(xp_design_2E$batch3)
write.table(xp_design_2E, "xp_design_2E.csv", sep=",", row.names=TRUE)

#make DESeq dataset
suppressPackageStartupMessages(library(DESeq2))
counts_E <- data.matrix(counts_E, rownames.force=T)
class(counts_E)

xp_design_2E$group <- rep(c("e02","e06","e24","m02","m06","m24","n02","n06","n24"), c(7,5,5,5,5,5,5,5,5))
dds <- DESeqDataSetFromMatrix(countData = counts_E,
                              colData = xp_design_2E,
                              design = ~ batch1 + batch2 + batch3 + time + treatment
)
dds$treatment <- relevel(dds$treatment, ref = "e")

#normalise and transform data
library(ggplot2)
library(matrixStats)
dds_norm <- estimateSizeFactors(dds)
counts_normalised <- counts(dds_norm, normalized = TRUE)

#create dds_norm_vst 
dds_vst <- vst(dds)
dds_norm_vst <- vst(dds_norm)
counts_norm_vst <- assay(dds_norm_vst, normalized = TRUE)
counts_vst <- assay(dds_vst, normalized = FALSE)

#PCA using the plotPCA function variance of treatment, shows there is no batch effect and E02_1_1 and E02_1_2 are close together
plotPCA(dds_vst, intgroup = c("treatment"))
z <- plotPCA(dds_vst, intgroup = c("time", "treatment"), returnData = TRUE)
ggplot(data = z, aes(PC1, PC2, colour = treatment)) +
  geom_point(size=3) +
  geom_text(aes(label=name), size=3, position = position_nudge(x=0.7)) +
  scale_fill_manual(values = clr,
                    guide = 'none')
z <- plotPCA(dds_vst, intgroup = c("batch3", "treatment"), returnData = TRUE)
ggplot(data = z, aes(PC1, PC2, colour = batch3)) +
  geom_point(size=3) +
  geom_text(aes(label=name), size=3, position = position_nudge(x=0.7)) +
  scale_fill_manual(values = clr,
                    guide = 'none')
z <- plotPCA(dds_vst, intgroup = c("batch2", "treatment"), returnData = TRUE)
ggplot(data = z, aes(PC1, PC2, colour = batch2)) +
  geom_point(size=3) +
  geom_text(aes(label=name), size=3, position = position_nudge(x=0.7)) +
  scale_fill_manual(values = clr,
                    guide = 'none')
z <- plotPCA(dds_vst, intgroup = c("batch1", "treatment"), returnData = TRUE)
ggplot(data = z, aes(PC1, PC2, colour = batch1)) +
  geom_point(size=3) +
  geom_text(aes(label=name), size=3, position = position_nudge(x=0.7)) +
  scale_fill_manual(values = clr,
                    guide = 'none')

