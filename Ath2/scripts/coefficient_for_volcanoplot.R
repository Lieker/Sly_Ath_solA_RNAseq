library(DESeq2)
source("scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "inputs/counts.csv",
                        xp_design_csv_file = "inputs/xp_design.csv",
                        trtm = c("a","b"),
                        ref_treatment = "a",
                        treatment2 = "b") {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2)
  r <- resultsNames(dds)[2]
  return(r)
}