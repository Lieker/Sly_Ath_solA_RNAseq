library(vegan)
source("scripts/filter_counts_vst_based_on_organ.R")
source("scripts/get_xp_design.R")

do_permanova <- function(xp_design_csv_file = "inputs/xp_design.csv",
                         counts_csv_file = "inputs/raw_counts.csv",
                         plantpart = "root") {
  
  xp_design <- get_xp_design(xp_design_csv_file)
  xp_design_plantpart <- xp_design[xp_design$organ == plantpart,]
  
  filtered_counts_vst_t <- filter_counts_vst_based_on_organ(counts_csv_file, plantpart)
  
  x <- adonis(filtered_counts_vst_t~treatment, data=xp_design_plantpart, permutations = 9999, method = "euclidean")
  return(x)
}
