suppressPackageStartupMessages(library(tidyverse))
source("scripts/pca_function.R")

produce_score_df <- function(counts = filtered_counts, 
                             xp_design_csv_file = "xp_design.csv",... 
#                             pc_x_axis. = pc_x_axis, 
#                             pc_y_axis. = pc_y_axis)  
 {
  
  pca_results <- mypca(x = counts, center = TRUE, scale = TRUE)
  scores <- pca_results$score
  
  # keep only components of interest (numbers here to filter the score df)
  first_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis.))
  second_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis.))
  
  scores_filtered <- scores[,c(first_pc, second_pc)]
    
  # add xp_design extra info for future plot
  xp_design = read.csv(file = xp_design_csv_file, header = TRUE) 
  
  scores_with_xp_design_info <- scores_filtered %>% 
    as.data.frame() %>% 
    rownames_to_column("sample") %>% 
    left_join(x = ., y = xp_design, by = "sample") %>% 
    column_to_rownames("sample")  
  
  return(scores_with_xp_design_info)
}