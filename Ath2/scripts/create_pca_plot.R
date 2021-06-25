plot_pca <- function(count_csv_file = "Ath2/input/counts.csv",
                     xp_design_csv_file = "Ath2/input/xp_design.csv",
                     pc_x_axis = "PC1", 
                     pc_y_axis = "PC2",
                     trm = c("a","b","c","d","e","f","g"),
                     pca_colour = "treatment",
                     pca_shape = 'solA') {
  
  source("Ath2/scripts/produce_scaled_counts_matrix.R")
  source("Ath2/scripts/filter_counts_based_on_treatment.R")
  source("Ath2/scripts/extract_variance.R")
  
  produce_score_df <- function(counts = filtered_counts, 
                               xp_design_csv_file) {
    
    pca_results <- mypca(x = counts, center = TRUE, scale = TRUE)
    scores <- pca_results$score
    
    # keep only components of interest (numbers here to filter the score df)
    first_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis))
    second_pc <- as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis))
    
    scores_filtered <- scores[,c(first_pc, second_pc)]
    
    # add xp_design extra info for future plot
    xp_design = read.csv(file = xp_design_csv_file, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM") 
    
    scores_with_xp_design_info <- scores_filtered %>% 
      as.data.frame() %>% 
      rownames_to_column("sample") %>% 
      left_join(x = ., y = xp_design, by = "sample") %>% 
      column_to_rownames("sample")  
    
    return(scores_with_xp_design_info)
  }
  
  
  
  # scale raw counts the DESeq2 way
  scaled_counts = produce_scaled_counts_matrix(count_csv_file = count_csv_file,
                                               xp_design_csv_file = xp_design_csv_file)
  
  # based on a timepoint, filter the corresponding scaled_counts matrix
  filtered_counts <- filter_counts_based_on_treatment(counts = scaled_counts, 
                                                      xp_design_csv_file = xp_design_csv_file, 
                                                      tr = trm) 
  
  # compute PCA and return scores as a dataframe with additional XP info
  # also returns explained variance per component 
   score_df <- produce_score_df(counts = filtered_counts, 
                               xp_design_csv_file = xp_design_csv_file)
  
  
  explained_variance_per_component <- extract_explained_variance_per_component(counts = filtered_counts) 
  
  # PCA plot
  variance_pc_x <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_x_axis))]
  variance_pc_y <- explained_variance_per_component[as.integer(gsub(pattern = "PC",replacement = "", x = pc_y_axis))]
  
  score_df[,pca_shape] <- as.factor(score_df[,pca_shape])
 
  pca_plot <- ggplot(data = score_df) +
    geom_point(aes_string(x = pc_x_axis, 
                          y = pc_y_axis, 
                          col = pca_colour,
                          shape = pca_shape, 
                          size = 3)) +
    theme_minimal() +
    scale_color_manual(values=c('#0099FF','#E69F00', '#56B4E9',"#e98b56","#56e9d5","#566be9","#e956b4")) +
    guides(size = FALSE) +
    xlab(paste0(pc_x_axis," ", variance_pc_x, "%")) +
    ylab(paste0(pc_y_axis," ", variance_pc_y, "%"))+
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15))

  return(pca_plot)
}


