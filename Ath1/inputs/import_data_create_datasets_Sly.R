library(dplyr)
library(stringr)

setwd("C:/Users/levla/github/Ath_RNAseq/inputs")


#-----------------------------------------------------------------------------------------------------

#make dds dataset for complete dataset (all samples included)
counts_raw <- read.delim("counts_Sly.txt", header = F, stringsAsFactors = F)
counts <- counts_raw[-c(1),]
headers <- counts[1,]
counts <- counts[-c(1),]
genes <- counts[,1]
counts <- counts[,-1]
genes <- str_replace(genes, "gene:","")
row.names(counts) <- genes
headers <- str_replace(headers, "_Aligned.sortedByCoord.out.bam","")
headers <- str_replace(headers, "mapped/","")
headers <- headers[2:26]
headers <- str_replace(headers, "_R1_001_qc", "")
names(counts) <- headers
counts_raw_clean <- counts
counts <- counts[,-c(1,2,3,4,5)]

#compile the xp_design table
treatment <- rep(c("no_solA","yes_solA"), c(10,10))
organ <- rep(c("root","shoot","root","shoot"), c(5,5,5,5))
headers <- colnames(counts)
xp_design <- data.frame(headers, treatment, organ)
colnames(xp_design) <- c("sample","treatment","organ")

#check if counts and xp_design match
all(colnames(counts) %in% xp_design$sample)
all(colnames(counts) == xp_design$sample)
xp_design$treatment <- as.factor(xp_design$treatment)
xp_design$organ <- as.factor(xp_design$organ)
write.table(xp_design, "xp_design.csv", sep=",", row.names=TRUE)
write.table(counts, "raw_counts.csv", sep=",", row.names=TRUE)
