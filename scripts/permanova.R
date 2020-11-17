library(vegan)
source("scripts/produce_vst_norm_counts_matrix.R")
source("scripts/get_xp_design.R")

do_permanova <- function(counts_vst_t = counts_vst_t,
                         xp_des = xp_design) {
  adonis(counts_vst_t~treatment, data=xp_des, permutations = 9999, method = "euclidean")
}

