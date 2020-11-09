create_pca_plot <- function(counts = scaled_counts, 
                            design = xp_design, 
                            time = NULL, 
                            pc_x_axis = 1,
                            pc_y_axis = 2,
                            center = TRUE, 
                            scale = TRUE) {
  
  source("scripts/pca_function.R")
  
  # no time filter. PCA on the whole set
  if (is.null(time)){
    pca_res <- mypca(counts, center = center, scale = scale)

    # Select components of interest
    scores = pca_res$scores[,c(pc_x_axis,pc_y_axis)]
    explained_variance = pca_res$explained_variance[c(pc_x_axis,pc_y_axis),]
  
    # add experimental design info (samples to factors)
    head(scores)
    
    # plot PCA
    # plot the scores of the selected 2 components
    
    
    p <- ggplot(scores) + 
      geom_point(aes(x = paste0("PC",pc_x_axis), 
                     y = paste0("PC",pc_y_axis), 
                     shape = time)) + 
      xlab(paste0('PC1(',explained_variance[1],'%)')) + 
      ylab(paste0('PC2(',explained_variance[2],'%)')) + 
      ggtitle('PCA score plot')
    p
    
  }
}
  