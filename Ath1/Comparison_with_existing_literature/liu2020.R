alldegs <- read.csv("alldegs.csv")
liu <- read.csv("Scheible_data/liu_DEGs.csv")
names(liu) <- c("genes", "lfcliu2020","qvalliu2020")
mergeliu <- merge(alldegs, liu, by = "genes")
mergeliu <- mergeliu[,-c(2,4,5,6,9,10,13,16)]
write.csv(mergeliu, "Scheible_data/DEG_incommonwith_Liu2020.csv", row.names = F)
