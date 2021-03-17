library(mixOmics)
source("scripts/filter_counts_vst_based_on_treatment.R")
source("scripts/get_xp_design.R")

do_plsda_perf <- function(xp_design_csv_file = "inputs/xp_design.csv",
                     counts_csv_file = "inputs/counts.csv",
                     trtm = c("a","b","c","d","e","f","g"),
                     n = 999) 
  {
  
  xp_design <- get_xp_design(xp_design_csv_file = xp_design_csv_file)
  xp_design <- xp_design[which(xp_design$treatment %in% trtm),]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_treatment(counts_csv_file, trtm)
  
  p <- plsda(filtered_counts_vst_t, xp_design$treatment, ncomp = 10)
  
  perf.plsda <- perf(p, 
                     validation = "Mfold", 
                     folds = 5, 
                     progressBar = FALSE, 
                     dist = 'max.dist', 
                     auc = FALSE, 
                     nrepeat = n)

  nc <- perf.plsda$choice.ncomp
  return(nc) #select the 'overall' ideal number of ideal components (instead of BER)
}
