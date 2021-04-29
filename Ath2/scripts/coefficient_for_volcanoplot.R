library(DESeq2)
source("Ath2/scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "Ath2/input/counts.csv",
                        xp_design_csv_file = "Ath2/input/xp_design.csv",
                        trtm = c("a","b"),
                        ref_treatment = "a",
                        treatment2 = "b",
                        method = "treatment" #this parameter chooses which formula design will be chosen: ~treatment or ~solA
                        ) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       trtm,
                       ref_treatment,
                       treatment2,
                       method)
  if(method == "treatment"){
    r <- resultsNames(dds)[2]
  } else if(method == "solA") {
    r <- resultsNames(dds)[4]
    dds$condition <- relevel(dds$solA, ref = "no")
  } else if(method == "N:solA+solA") {
    r <- resultsNames(dds)[5]
    dds$condition <- relevel(dds$solA, ref = "no")
  } else if(method == "N:solA"){
    r <- resultsNames(dds)[5]
    dds$condition <- relevel(dds$solA, ref = "no")
  }
    
  
  return(r)
}
