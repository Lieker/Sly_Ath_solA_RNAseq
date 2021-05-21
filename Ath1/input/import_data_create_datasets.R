library(dplyr)
library(stringr)

setwd("Z:/Experiments/Arabidopsis_RNAseq/R")


#-----------------------------------------------------------------------------------------------------

#make dds dataset for complete dataset (all samples included)
counts_raw <- read.csv("countsc.csv", header = T, stringsAsFactors = F)
counts <- counts_raw[-c(1),]
headers <- counts[1,]
counts <- counts[-c(1),]
genes <- counts[,1]
counts <- counts[,-1]
row.names(counts) <- genes
headers <- str_replace(headers, "_Aligned.sortedByCoord.out.bam","")
headers <- str_replace(headers, "mapped/","")
headers <- headers[2:52]
names(counts) <- headers
counts_raw_clean <- counts
counts <- counts[,-c(1,2,3,4,5)]

#compile the xp_design table
counts <- counts %>% select(sort(current_vars()))
batches <- read.csv("batches.csv", header = T, stringsAsFactors = F)
treatment <- rep(c("e","m","n"), c(16, 15, 15))
time <- rep(c("2","6","24","2","6","24","2","6","24"), c(6,5,5,5,5,5,5,5,5))
headers <- colnames(counts)
xp_design <- data.frame(headers, treatment, time, batches$batch, batches$batch2, batches$batch3)
colnames(xp_design) <- c("sample","treatment","time","batch1","batch2","batch3")

#check if counts and xp_design match
all(colnames(counts) %in% xp_design$sample)
all(colnames(counts) == xp_design$sample)
xp_design$batch1 <- as.factor(xp_design$batch1)
xp_design$treatment <- as.factor(xp_design$treatment)
xp_design$batch3 <- as.factor(xp_design$batch3)
write.table(xp_design, "xp_design.csv", sep=",", row.names=TRUE)

#make dds dataset
suppressPackageStartupMessages(library(DESeq2))
counts <- mutate_all(counts, function(x) as.numeric(as.character(x)))
row.names(counts) <- genes
counts <- data.matrix(counts, rownames.force=T)
class(counts)
counts[1,1] == counts_raw[3,7]

xp_design$group <- rep(c("e02","e06","e24","m02","m06","m24","n02","n06","n24"), c(6,5,5,5,5,5,5,5,5))
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = xp_design,
                              design = ~ batch1 + batch2 + batch3 + time + treatment
)
dds$treatment <- relevel(dds$treatment, ref = "e")

#-------------------------------------------------------------------------------------------

#make datasets for all individual timepoints

#make dds for t=2, t=6 and t=24
xp_design_2 <- xp_design %>% filter(time == 2)
xp_design_6 <- xp_design %>% filter(time == 6)
xp_design_24 <- xp_design %>% filter(time == 24)
counts_2 <- counts[,colnames(counts) %in% xp_design_2$sample]
counts_6 <- counts[,colnames(counts) %in% xp_design_6$sample]
counts_24 <- counts[,colnames(counts) %in% xp_design_24$sample]
dds_2 <- DESeqDataSetFromMatrix(countData = counts_2,
                                colData = xp_design_2,
                                design = ~batch1 + batch2 + batch3 + treatment)
dds_6 <- DESeqDataSetFromMatrix(countData = counts_6,
                                colData = xp_design_6,
                                design = ~ batch2 + treatment)
dds_24 <- DESeqDataSetFromMatrix(countData = counts_24,
                                 colData = xp_design_24,
                                 design = ~batch2 + treatment)
dds_2$treatment <- relevel(dds_2$treatment, ref = "e")
dds_6$treatment <- relevel(dds_6$treatment, ref = "e")
dds_24$treatment <- relevel(dds_24$treatment, ref = "e")


#-------------------------------------------------------------------------------------------

#make datasets for individual timepoints, M and N treatments separated

#create dds only for 2m vs e
xp_design_2m <- xp_design %>% filter(time == 2 & treatment != "n")
counts_2m <- counts[,colnames(counts) %in% xp_design_2m$sample]

dds_2m <- DESeqDataSetFromMatrix(countData = counts_2m,
                                 colData = xp_design_2m,
                                 design = ~batch1 + batch2 + treatment)
dds_2m$treatment <- relevel(dds_2m$treatment, ref = "e")

#create dds only for 2n vs e
xp_design_2n <- xp_design %>% filter(time == 2 & treatment != "m")
counts_2n <- counts[,colnames(counts) %in% xp_design_2n$sample]

dds_2n <- DESeqDataSetFromMatrix(countData = counts_2n,
                                 colData = xp_design_2n,
                                 design = ~batch1 + treatment)
dds_2n$treatment <- relevel(dds_2n$treatment, ref = "e")

#create dds only for 6m vs e
xp_design_6m <- xp_design %>% filter(time == 6 & treatment != "n")
counts_6m <- counts[,colnames(counts) %in% xp_design_6m$sample]

dds_6m <- DESeqDataSetFromMatrix(countData = counts_6m,
                                 colData = xp_design_6m,
                                 design = ~batch2 + treatment)
dds_6m$treatment <- relevel(dds_6m$treatment, ref = "e")

#create dds only for 6n vs e
xp_design_6n <- xp_design %>% filter(time == 6 & treatment != "m")
counts_6n <- counts[,colnames(counts) %in% xp_design_6n$sample]

dds_6n <- DESeqDataSetFromMatrix(countData = counts_6n,
                                 colData = xp_design_6n,
                                 design = ~batch2 + treatment)
dds_6n$treatment <- relevel(dds_6n$treatment, ref = "e")

#create dds only for 24m vs e
xp_design_24m <- xp_design %>% filter(time == 24 & treatment != "n")
counts_24m <- counts[,colnames(counts) %in% xp_design_24m$sample]

dds_24m <- DESeqDataSetFromMatrix(countData = counts_24m,
                                  colData = xp_design_24m,
                                  design = ~batch2 + treatment)
dds_24m$treatment <- relevel(dds_24m$treatment, ref = "e")

#create dds only for 24n vs e
xp_design_24n <- xp_design %>% filter(time == 24 & treatment != "m")
counts_24n <- counts[,colnames(counts) %in% xp_design_24n$sample]

dds_24n <- DESeqDataSetFromMatrix(countData = counts_24n,
                                  colData = xp_design_24n,
                                  design = ~batch2 + treatment)
dds_24n$treatment <- relevel(dds_24n$treatment, ref = "e")