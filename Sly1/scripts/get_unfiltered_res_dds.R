library(DESeq2)


get_unfiltered_res_dds <- function(counts_csv_file = "Sly1/input/counts.csv",
                                   xp_design_csv_file = "Sly1/input/xp_design.csv",
                                   plantpart = "root",
                                   ref_treatment = "no_solA",
                                   treatment2 = "yes_solA",
                                   method) {
  source("Sly1/scripts/get_DESeq_dds.R")
  
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       plantpart,
                       ref_treatment,
                       treatment2,
                       method)
  
  res <- results(dds)
  return(res)
}
