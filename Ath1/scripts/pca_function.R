mypca <- function(x, center = T, scale = T){  
  
  # remove columns containing only 0 values
  # not informative + cause svd() error
  x_without_zero_columns_matrix <- x[vapply(x,
                                     function(z) length(unique(z)) > 1,
                                     logical(1L))] %>% data.matrix(., rownames.force = T)
  x_without_zero_columns_frame <- x[vapply(x,
                                           function(z) length(unique(z)) > 1,
                                           logical(1L))]
  median_all_genes <- median(x_without_zero_columns_matrix)
  median_per_gene <- apply(x_without_zero_columns_frame, 2, median)
  x_without_zero_columns_frame_t <- x_without_zero_columns_frame %>% t() %>% as.data.frame()
  x_without_zero_columns <- x_without_zero_columns_frame_t %>% filter(median_per_gene > median_all_genes) %>% t() %>% as.data.frame()
  
  # perform SVD
  SVD <- svd(scale(x_without_zero_columns, center = center, scale = scale))
  
  # create scores data frame
  scores <- as.data.frame(SVD$u %*% diag(SVD$d))
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste0("PC", c(1:dim(scores)[2]))
  
  # create loadings data frams
  loadings <- data.frame(SVD$v)
  colnames(loadings) <- paste0("PC", c(1:dim(loadings)[2]))
  row.names(loadings) <- colnames(x_without_zero_columns)
  
  # create data frame for explained variances
  explained_var <- as.data.frame(round((SVD$d^2) / sum(SVD$d^2)*100, digits = 1))
  rownames(explained_var) <- paste0("PC", c(1:dim(loadings)[2]))
  colnames(explained_var) <- "exp_var"
  explained_var = mutate(explained_var, PC = 1:nrow(explained_var))
  
  # return results
  
  return (list("scores" = scores, "loadings" = loadings, "explained_var" = explained_var))
}
