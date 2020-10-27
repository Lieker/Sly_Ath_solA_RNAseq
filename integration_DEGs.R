#-------------------------compare DEGs of all treatment*timepoints -------------------------------------------

#do all other DEG analyses of other timepoints and treatments in separate scripts
#THEN
#check if DEGs of timepoints are similar
diff_2mn <- merge(diff_2m_annotated, diff_2n_annotated, by="ensembl_gene_id")
diff_2mn <- res_2mn_f[,-c(2,4:6,8,10:12)]
diff_24mn <- merge(diff_24m_annotated, diff_24n_annotated, by="ensembl_gene_id")
diff_24mn <- res_24mn_f[,-c(2,4:6,8,10:12)]

#these lists combined do not give any shared genes between all timepoints
write.table(diff_24mn, file = "DEG_shared_24mn.csv", sep=";", dec=".", row.names = FALSE)

#combine genes present in all samples that have 'oxygen' reaction
diff_oxygen <- merge(diff_2n_annotated, diff_6n_annotated, by='ensembl_gene_id')
diff_oxygen <- merge(diff_oxygen, diff_2mn, by = 'ensembl_gene_id')
#result is only 2 genes

alldegs <- read.csv("all_DEGs_trimmed.csv", header = T, stringsAsFactors = F)
alldegscount <- alldegs %>%
  group_by(ensembl_gene_id) %>%
  mutate(count = n())
multdegs <- alldegscount[alldegscount$count > 1, ]
write.table(multdegs, "multdegs.csv", sep=";", dec=".", row.names = TRUE)
