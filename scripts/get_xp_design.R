get_xp_design <- function(xp_design_csv_file = "xp_design.csv") {
  xp_design <- read.csv(file = xp_design_csv_file, 
                        header = TRUE, 
                        stringsAsFactors = FALSE, 
                        check.names = FALSE, fileEncoding = "UTF-8-BOM")
  return(xp_design)
}