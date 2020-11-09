source("scripts/create_pca_plot.R")

plot_pca(count_csv_file = "raw_counts.csv", 
         xp_design_csv_file = "xp_design.csv", 
         timepoint = 2, 
         pc_x_axis = "PC1", 
         pc_y_axis = "PC2", 
         pca_colour = "sequencing_batch")
