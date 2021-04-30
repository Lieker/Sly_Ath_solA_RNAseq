library(dplyr)
library(stringr)

setwd("C:/Users/levla/github/Ath_RNAseq/Ath2")


#-----------------------------------------------------------------------------------------------------

#make dds dataset for complete dataset (all samples included)
counts_raw2 <- read.delim("Ath2/input/counts.txt", header = F, stringsAsFactors = F) %>% t() %>% as.data.frame() 
row.names(counts_raw2) <- NULL
counts_raw2 <- counts_raw2 %>% column_to_rownames("V1")
counts_raw2 <- counts_raw2 %>% t() %>% as.data.frame()
row.names(counts_raw2) <- NULL
Geneid <- counts_raw2$Geneid
names(counts_raw2)[names(counts_raw2) == "Geneid"] <- "Transcriptid"
genes <- substr(genes, start = 1, stop = 9)

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
