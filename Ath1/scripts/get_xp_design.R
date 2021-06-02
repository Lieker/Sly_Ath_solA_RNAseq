get_xp_design <- function(xp_design_csv_file = "Ath1/input/xp_design.csv") {
  xp_design <<- read.csv(file = xp_design_csv_file, 
                         header = TRUE, 
                         stringsAsFactors = FALSE, 
                         check.names = FALSE, 
                         fileEncoding = "UTF-8-BOM")
  xp_design$time <- as.factor(xp_design$time)
  return(xp_design)
}
