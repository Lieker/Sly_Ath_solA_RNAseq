library(DESeq2)
source("scripts/get_DESeq_dds.R")

get_unfiltered_res_dds <- function(counts_csv_file = "inputs/counts.csv",
                                   xp_design_csv_file = "inputs/xp_design.csv",
                                   trtm = c("a","b"),
                                   ref_treatment = "a",
                                   treatment2 = "b") {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2)
  res <- results(dds)
  return(res)
}