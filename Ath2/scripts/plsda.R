library(mixOmics)
source("scripts/filter_counts_vst_based_on_treatment.R")
source("scripts/get_xp_design.R")
source("scripts/perf_plsda.R")

do_plsda <- function(xp_design_csv_file = "inputs/xp_design.csv",
                     counts_csv_file = "inputs/counts.csv",
                     trtm = c("a","b","c","d","e","f","g"),
                     nr = 999) {
  
  xp_design_f <- xp_design[which(xp_design$treatment %in% trtm),]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_treatment(counts_csv_file, trtm)
  
  nc <- do_plsda_perf(xp_design_csv_file, counts_csv_file, trtm, n = nr)
  
  q <- plsda(filtered_counts_vst_t, xp_design_f$treatment, ncomp = nc)
  
  r <- plotIndiv(q,
                 comp = 1:2,
                 group = xp_design_f$treatment,
                 ind.names = FALSE,
                 legend = TRUE,
                 ellipse = TRUE)
  return(r)
}
