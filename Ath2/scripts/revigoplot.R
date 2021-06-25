plot_dot <- function(x = "Ath2/output/panther/GOann-pantherfdr_replete.csv",
                     sign = 0.05)
  {
  gos <- read.csv(x, stringsAsFactors = FALSE)
  names(gos) <- c("GOterm","ref","ngenes","expected","over/under","foldenrich","rawPvalue","fdr")
  gos <- gos %>% dplyr::filter(gos$fdr < sign)
  p <- ggplot(gos, 
              aes(x=reorder(GOterm, 
                            foldenrich), 
                  y=foldenrich, 
                  size=ngenes, 
                  color=fdr)) + 
    geom_point(alpha = 0.8) + 
    theme_minimal() + 
    coord_flip() +
    scale_color_gradient(low = "firebrick1",
                         high = "blue4",
                         space = "Lab", 
                         limit = c(0, sign)) +
    theme(axis.title.y=element_blank())
  return(p)
  
}
