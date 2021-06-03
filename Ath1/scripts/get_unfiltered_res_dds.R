library(DESeq2)
source("Ath1/scripts/get_DESeq_dds.R")

get_unfiltered_res_dds <- function(counts_csv_file = "Ath1/input/counts.csv",
                                   xp_design_csv_file = "Ath1/input/xp_design.csv",
                                   trtm = c("ethanol","millimolar_solanoeclepinA"),
                                   ref_treatment = "ethanol",
                                   treatment2 = "millimolar_solanoeclepinA",
                                   method = "time+treatment",
                                   tp = c(2, 6, 24)) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method,
                       tp)
  res <- results(dds)
  return(res)
}
