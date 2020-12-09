nguyen <- read.csv("nguyen.csv", header = TRUE)
names(nguyen) <- c("genes", "logFC_nguyen2017", "Padj_nguyen2017")
degs_n <- left_join(degs, nguyen, by = "genes")
write.csv(degs_n, "allDEGs_s_n.csv", row.names = FALSE)
