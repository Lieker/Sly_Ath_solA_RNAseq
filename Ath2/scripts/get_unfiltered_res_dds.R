library(DESeq2)
source("scripts/get_DESeq_dds")

get_unfiltered_res_dds <- function(counts_csv_file = "raw_counts.csv",
                                   xp_design_csv_file = "xp_design.csv",
                                   timepoint = 2,
                                   ref_treatment = "ethanol",
                                   treatment2 = "millimolar_solanoeclepinA") {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       timepoint,
                       ref_treatment,
                       treatment2)
  res <- results(dds)
  return(res)
}
r<-get_unfiltered_res_dds()