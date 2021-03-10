library(DESeq2)


get_unfiltered_res_dds <- function(counts_csv_file = "inputs/raw_counts.csv",
                                   xp_design_csv_file = "inputs/xp_design.csv",
                                   plantpart = "root",
                                   ref_treatment = "no_solA",
                                   treatment2 = "yes_solA") {
  source("scripts/get_DESeq_dds.R")
  
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       plantpart,
                       ref_treatment,
                       treatment2)
  
  res <- results(dds)
  return(res)
}
