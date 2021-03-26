library(mixOmics)
source("scripts/filter_counts_based_on_compartment.R")
source("scripts/get_xp_design.R")

do_plsda_perf <- function(xp_design_csv_file = "inputs/xp_design.csv",
                     counts_csv_file = "inputs/raw_counts.csv",
                     plantpart = "root",
                     n = 999) 
  {
  
  xp_design_f <- xp_design[xp_design$organ == plantpart,]
  
  filtered_counts_vst_t <- filter_counts_based_on_compartment(counts_csv_file, comprtmt = plantpart)
  
  p <- plsda(filtered_counts_vst_t, xp_design_f$treatment, ncomp = 5)
  
  perf.plsda <- perf(p, 
                     validation = "Mfold", 
                     folds = 5, 
                     progressBar = FALSE, 
                     dist = 'max.dist', 
                     auc = FALSE, 
                     nrepeat = n)

  nc <- perf.plsda$choice.ncomp
  return(nc[1,1]) #select the 'overall' ideal number of ideal components (instead of BER)
}
