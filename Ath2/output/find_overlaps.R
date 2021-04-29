library(dplyr)


setwd("C:/Users/levla/github/Ath_RNAseq/Ath2")
Ath2_ab <- read.csv("outputs/annotated_DEGslist_ab.csv", sep = ";", stringsAsFactors = F)
Ath2_cd <- read.csv("outputs/annotated_DEGslist_cd.csv", sep = ";", stringsAsFactors = F)
Ath2_ef <- read.csv("outputs/annotated_DEGslist_ef.csv", sep = ";", stringsAsFactors = F)
Ath2_overlap_abcd <- inner_join(Ath2_ab, Ath2_cd, by = "genes", suffix = c(".ab",".cd"), copy = F)
Ath_overlap_abcdef <- inner_join(Ath2_overlap_abcd, Ath2_ef, by = "genes", suffix = c(".abcd", ".ef"), copy = F)
write.csv(Ath_overlap_abcdef, "outputs/overlap_degs_within_Ath2.csv", row.names = F)
setwd("C:/Users/levla/github/Ath_RNAseq/Sly1")
Slyr <- read.csv("outputs/annotated_DEGslist_root.csv", sep = ";", stringsAsFactors = F)
Slys <- read.csv("outputs/annotated_DEGslist_shoot.csv", sep = ";", stringsAsFactors = F)
names(Ath_overlap_abcdef)[1] <- "athaliana_eg_homolog_ensembl_gene"
all2 <- inner_join(Ath_overlap_abcdef, Slyr, by = "athaliana_eg_homolog_ensembl_gene", suffix = c("Ath", "Sl"))
all2 <- inner_join(Ath_overlap_abcdef, Slys, by = "athaliana_eg_homolog_ensembl_gene", suffix = c("Ath", "Sl"))
names(Ath2_ab)[1] <- "athaliana_eg_homolog_ensembl_gene"
names(Ath2_cd)[1] <- "athaliana_eg_homolog_ensembl_gene"
names(Ath2_ef)[1] <- "athaliana_eg_homolog_ensembl_gene"
some2 <- inner_join(Ath2_cd, Slyr, by = "athaliana_eg_homolog_ensembl_gene", suffix = c("A","Sr"))
Ath1 <- read.csv("outputs/Ath1_alldegslist.csv", sep = ";", stringsAsFactors = F)
names(Ath1)[1] <- "athaliana_eg_homolog_ensembl_gene"
ath1sly <- inner_join(Ath1, Slyr, by = "athaliana_eg_homolog_ensembl_gene", suffix = c("A1", "Sr"))
