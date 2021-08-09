plotcounts <- function(gene, d = dds, xpd = xp_design, trtm = c("replete","Pstarv","NPstarv"),
                       name = gene) {
  data <- plotCounts(d, gene,
                     intgroup = "solA", 
                     returnData = TRUE) %>% rownames_to_column("sample") 
  data <- left_join(data, xpd, by = c("sample", "solA"))
  tt <- c(rep("replete",10), rep("Pstarv",10), rep("NPstarv",10),rep("Nstarv",5))
  
  data$treatment <- tt
  data$treatment <- factor(data$treatment, levels = c("replete","Pstarv","NPstarv","Nstarv"))
  data <- data[order(data$treatment),]
  data <- data %>% dplyr::filter(treatment %in% trtm)
  
  g <- ggplot(data, aes(x=treatment, y=count)) + 
    geom_jitter(aes(color = solA, shape = solA), 
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0), size = 2) +
    theme_classic() +
    ggtitle(paste0(name)) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
    scale_color_manual(values = c("#E41A1C","#377EB8"), name = "Treatment", labels = c("control", "solA")) +
    scale_shape_discrete(name = "Treatment", labels = c("control", "solA")) +
    xlab("")
  return(g)
}
