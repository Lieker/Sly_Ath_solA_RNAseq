source("scripts/produce_scaled_counts_matrix.R")
source("scripts/filter_counts_based_on_time.R")
source("scripts/produce_score_df.R")
source("scripts/extract_variance.R")

plot_pca <- function(count_csv_file = "raw_counts.csv",
                     xp_design_csv_file = "xp_design.csv",
                     pc_x_axis = "PC1", 
                     pc_y_axis = "PC2",
                     timepoint = 2,
                     pca_colour = "rna_isolation_batch") {
  
  # scale raw counts the DESeq2 way
  scaled_counts = produce_scaled_counts_matrix(count_csv_file = count_csv_file,
                                               xp_design_csv_file = xp_design_csv_file)
  
  # based on a timepoint, filter the corresponding scaled_counts matrix
  filtered_counts <- filter_counts_based_on_time(counts = scaled_counts, 
                                                 xp_design_csv_file = xp_design_csv_file, 
                                                 time = timepoint) 
  
  # compute PCA and return scores as a dataframe with additional XP info
  # also returns explained variance per component 
  score_df <- produce_score_df(counts = filtered_counts, 
                               xp_design_csv_file = "xp_design.csv")
  
  
  explained_variance_per_component <- extract_explained_variance_per_component(counts = filtered_counts) 
  
  # PCA plot
  variance_pc_x <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis))]
  variance_pc_y <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis))]
  
  
  pca_plot <- ggplot(data = score_df) +
    geom_point(aes_string(x = pc_x_axis, 
                          y = pc_y_axis, 
                          col = pca_colour), 
               size = 3) +
    xlab(paste0(pc_x_axis," ", variance_pc_x, "%")) +
    ylab(paste0(pc_y_axis," ", variance_pc_y, "%"))

  return(pca_plot)
}


plot_pca(timepoint = 2)

