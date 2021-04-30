library(DESeq2)
source("Sly1/scripts/get_DESeq_dds.R")

coefficient <- function(counts_csv_file = "Sly1/input/raw_counts.csv",
                        xp_design_csv_file = "Sly1/input/xp_design.csv",
                        plantpart = "root",
                        ref_treatment = "no_solA",
                        treatment2 = "yes_solA",
                        method) {
  dds <- get_DESeq_dds(counts_csv_file,
                       xp_design_csv_file,
                       plantpart,
                       ref_treatment,
                       treatment2,
                       method)
  if(method == "treatment"){
    r <- resultsNames(dds)[2]
  } else if(method == "compartment+treatment"){
    r <- resultsNames(dds)[3]
  }
  
  return(r)
}