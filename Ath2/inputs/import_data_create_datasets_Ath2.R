library(dplyr)
library(stringr)

setwd("C:/Users/levla/github/Ath_RNAseq/Ath2")


#-----------------------------------------------------------------------------------------------------

#make dds dataset for complete dataset (all samples included)
counts_raw <- read.delim("inputs/counts_Ath2.txt", header = F, stringsAsFactors = F)
counts <- counts_raw[-c(1),]
headers <- counts[1,]
counts <- counts[-c(1),]
genes <- counts[,1]
counts <- counts[,-1]
genes <- str_replace(genes, "gene:","")
row.names(counts) <- genes
headers <- str_replace(headers, "_Aligned.sortedByCoord.out.bam","")
headers <- str_replace(headers, "mapped/","")
headers <- headers[7:41]
headers <- str_replace(headers, "_R1_001_qc", "")
counts_raw_clean <- counts
counts <- counts[,-c(1,2,3,4,5)]
names(counts) <- headers


#compile the xp_design table
treatment <- rep(c("a","b","c","d","e","f","g"), c(5,5,5,5,5,5,5))
solA <- rep(c("no","yes","no","yes","no","yes","no"), c(5,5,5,5,5,5,5))
N <- rep(c("yes","no"), c(20,15))
P <- rep(c("yes","no","yes"), c(10,20,5))
headers <- colnames(counts)
xp_design <- data.frame(headers, treatment, solA, N, P)
colnames(xp_design) <- c("sample","treatment","solA", "N", "P")

#check if counts and xp_design match
all(colnames(counts) %in% xp_design$sample)
all(colnames(counts) == xp_design$sample)
xp_design$treatment <- as.factor(xp_design$treatment)
xp_design$solA <- as.factor(xp_design$solA)
xp_design$N <- as.factor(xp_design$N)
xp_design$P <- as.factor(xp_design$P)
write.table(xp_design, "inputs/xp_design.csv", sep=",", row.names=F)
counts <- counts %>% rownames_to_column("Geneid")
write.table(counts, "inputs/counts.csv", sep=",", row.names=F)
