library(DESeq2)
source("scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "inputs/raw_counts.csv",
                        xp_design_csv_file = "inputs/xp_design.csv",
                        plantpart = "root",
                        ref_treatment = "no_solA",
                        treatment2 = "yes_solA") {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       plantpart,
                       ref_treatment,
                       treatment2)
  r <- resultsNames(dds)[2]
  return(r)
}