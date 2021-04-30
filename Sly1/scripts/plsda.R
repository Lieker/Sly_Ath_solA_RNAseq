source("Sly1/scripts/filter_counts_based_on_compartment.R")
source("Sly1/scripts/get_xp_design.R")
source("Sly1/scripts/perf_plsda.R")

do_plsda <- function(xp_design_csv_file = "Sly1/input/xp_design.csv",
                     counts_csv_file = "Sly1/input/raw_counts.csv",
                     plantpart = "root",
                     nc, #evaluate the performance of this pls-da to determine the nc
                     ttl = "") {
  xp_design <- get_xp_design()
  
  xp_design_f <- xp_design[xp_design$compartment %in% plantpart,]
  
  filtered_counts_vst_t <- filter_counts_based_on_compartment(comprtmt = plantpart) %>% data.matrix()
  
  q <- plsda(filtered_counts_vst_t, xp_design_f$treatment, ncomp = nc)
  
  r <- plotIndiv(q,
                 comp = 1:2,
                 group = xp_design_f$treatment,
                 ind.names = FALSE,
                 legend = TRUE,
                 ellipse = TRUE,
                 title = ttl)
  return(r)
}
