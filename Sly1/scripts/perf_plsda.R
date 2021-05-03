source("Sly1/scripts/filter_counts_based_on_compartment.R")
source("Sly1/scripts/get_xp_design.R")

do_plsda_perf <- function(xp_design_csv_file = "Sly1/input/xp_design.csv",
                     counts_csv_file = "Sly1/input/counts.csv",
                     plantpart = "root",
                     n = 9) 
  {
  xp_design <- get_xp_design(xp_design_csv_file)
  xp_design_f <- xp_design[xp_design$compartment %in% plantpart,]
  
  filtered_counts_vst_t <- filter_counts_based_on_compartment(comprtmt = plantpart) %>% data.matrix()

  p <- plsda(filtered_counts_vst_t, xp_design_f$treatment, ncomp = 5)
  
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
