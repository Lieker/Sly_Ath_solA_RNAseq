library(mixOmics)
source("scripts/filter_counts_vst_based_on_time.R")
source("scripts/get_xp_design.R")
source("scripts/perf_plsda.R")

do_plsda <- function(xp_design_csv_file = "inputs/xp_design.csv",
                     counts_csv_file = "inputs/raw_counts.csv",
                     timepoint = 2,
                     nr = 999) {
  
  xp_design_f <- xp_design[xp_design$time == timepoint,]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_time(counts_csv_file, timepoint)
  
  nc <- do_plsda_perf(xp_design_csv_file, counts_csv_file, timepoint, n = nr)
  
  q <- plsda(filtered_counts_vst_t, xp_design_f$treatment, ncomp = nc)
  
  background <- background.predict(q, comp.predicted=2, dist = "max.dist") 
  
  r <- plotIndiv(q,
                 comp = 1:2,
                 group = xp_design_f$treatment,
                 ind.names = FALSE,
                 legend = TRUE,
                 ellipse = TRUE,
                 background = background)
  return(r)
}

r <- do_plsda(nr = 9)
