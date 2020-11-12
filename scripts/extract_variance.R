source("scripts/pca_function.R")

extract_explained_variance_per_component <- function(counts = filtered_counts) {
  pca_results <- mypca(x = counts, center = TRUE, scale = TRUE)
  explained_variance <- pca_results$explained_var$exp_var
  
  return(explained_variance)
}