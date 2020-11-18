library(mixOmics)
source("scripts/filter_counts_vst_based_on_time.R")
source("scripts/get_xp_design.R")

do_plsda_perf <- function(xp_design_csv_file = "xp_design.csv",
                     counts_csv_file = "raw_counts.csv",
                     timepoint = 2,
                     n = 999) 
  {
  
  xp_design_f <- xp_design[xp_design$time == timepoint,]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_time(counts_csv_file, timepoint)
  
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