library(mixOmics)
source("scripts/filter_counts_based_on_compartment.R")
source("scripts/get_xp_design.R")
source("scripts/perf_plsda.R")

do_plsda <- function(xp_design_csv_file = "inputs/xp_design.csv",
                     counts_csv_file = "inputs/raw_counts.csv",
                     plantpart = "root",
                     nr = 999) {
  xp_design <- get_xp_design()
  
  xp_design_f <- xp_design[xp_design$organ == plantpart,]
  
  filtered_counts_vst_t <- filter_counts_based_on_compartment(counts_csv_file, comprtmt = plantpart) %>% data.matrix()
  
  nc <- do_plsda_perf(xp_design_csv_file, counts_csv_file, plantpart, n = nr)
  
  q <- plsda(filtered_counts_vst_t, xp_design_f$organ, ncomp = nc)
  
  r <- plotIndiv(q,
                 comp = 1:2,
                 group = xp_design_f$treatment,
                 ind.names = FALSE,
                 legend = TRUE,
                 ellipse = TRUE)
  return(r)
}
