plot_scree <- function(count_csv_file = "Sly1/input/counts.csv",
                       xp_design_csv_file = "Sly1/input/xp_design.csv",
                       comp = "root"
)
  {
  source("Sly1/scripts/filter_counts_based_on_compartment.R")
  source("Sly1/scripts/pca_function.R")
  scaled_counts = produce_scaled_counts_matrix(count_csv_file = count_csv_file,
                                               xp_design_csv_file = xp_design_csv_file)
    
  filtered_counts <- filter_counts_based_on_compartment(counts = scaled_counts, 
                                                      xp_design_csv_file = xp_design_csv_file, 
                                                      comprtmt = comp)
  pca_results <- mypca(x = filtered_counts)
  
  s <- ggplot(pca_results$explained_var, 
       aes(x = seq(from = 1, to = nrow(pca_results$explained_var)), 
           y = exp_var)) +
  ylab('explained variance (%)') + 
  ggtitle('Explained variance per component') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component number") +
  scale_x_continuous(breaks = seq(
    from = 1, 
    to = nrow(pca_results$explained_var)))
  
  return(s)
}
