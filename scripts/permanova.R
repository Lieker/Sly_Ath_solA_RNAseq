library(vegan)
source("scripts/filter_counts_vst_based_on_time.R")
source("scripts/get_xp_design.R")

do_permanova <- function(xp_design_csv_file = "xp_design.csv",
                         counts_csv_file = "raw_counts.csv",
                         timepoint = 2) {
  
  xp_design <- get_xp_design()
  xp_design_timepoint <- xp_design[xp_design$time == timepoint,]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_time(counts_csv_file, timepoint)
  
  x <- adonis(filtered_counts_vst_t~treatment, data=xp_design_timepoint, permutations = 999, method = "euclidean")
  return(x)
}

p <- do_permanova()
