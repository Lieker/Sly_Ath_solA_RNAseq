library(DESeq2)
source("scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "raw_counts.csv",
                        xp_design_csv_file = "xp_design.csv",
                        timepoint = 2,
                        ref_treatment = "ethanol",
                        treatment2 = "millimolar_solanoeclepinA") {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       timepoint,
                       ref_treatment,
                       treatment2)
  r <- resultsNames(dds)[2]
  return(r)
}