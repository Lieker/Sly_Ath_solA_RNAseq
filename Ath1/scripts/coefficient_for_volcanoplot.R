library(DESeq2)
source("Ath1/scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "Ath1/input/counts.csv",
                        xp_design_csv_file = "Ath1/input/xp_design.csv",
                        trtm = c("ethanol","millimolar_solanoeclepinA"),
                        ref_treatment = "ethanol",
                        treatment2 = "millimolar_solanoeclepinA",
                        method = "time+treatment", #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                        tp = c(2,6,24)) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method,
                       tp)
  if(method == "time+treatment"){
    r <- resultsNames(dds)[11]
    dds$condition <- relevel(dds$treatment, ref = "ethanol")
  } else if(method == "treatment"){
    r <- resultsNames(dds)[9]
    dds$condition <- relevel(dds$treatment, ref = "ethanol")
  }
  
  return(r)
}
