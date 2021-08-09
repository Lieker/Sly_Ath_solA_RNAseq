plotcounts <- function(gene, d = dds, genename,rs = "rs") {
  data <- plotCounts(d, gene,
                     intgroup = "group", 
                     returnData = TRUE) %>% rownames_to_column("sample") %>% separate(., 
                                                                                      col = group, 
                                                                                      into = c("compartment", "treatment"), 
                                                                                      sep = "t_", )
  data$e <- c(rep("t", nrow(data)))
  data$comp <- paste0(data$compartment,data$e)
  data$treatment <- factor(data$treatment, levels = c("no_solA", "yes_solA"), ordered = TRUE)
  data <- data[order(data$treatment), ]
  
  g <- ggplot(data, aes(x=comp, y=count)) + 
    geom_jitter(aes(shape = treatment, color = treatment), 
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0), size = 2) +
    theme_classic() +
    ggtitle(paste0(genename)) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
    theme(axis.title.x = element_blank()) +
    scale_color_manual(values = c("#E41A1C","#377EB8"), name = "Treatment", labels = c("control", "solA")) +
    scale_shape_discrete(name = "Treatment", labels = c("control", "solA"))
  
  if(rs == "r") {
    data <- data[data$comp == "root",]
    g <- ggplot(data, aes(x=comp, y=count)) + 
      geom_jitter(aes(shape = treatment, color = treatment), 
                  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0), size = 2) +
      theme_classic() +
      ggtitle(paste0(genename)) +
      scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
      theme(axis.title.x = element_blank()) +
      scale_color_manual(values = c("#E41A1C","#377EB8"), name = "Treatment", labels = c("control", "solA")) +
      scale_shape_discrete(name = "Treatment", labels = c("control", "solA"))
  } else if(rs == "s") {
    data <- data[data$comp == "shoot",]
    g <- ggplot(data, aes(x=comp, y=count)) + 
      geom_jitter(aes(shape = treatment, color = treatment), 
                  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0), size = 2) +
      theme_classic() +
      ggtitle(paste0(genename)) +
      scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000)) +
      theme(axis.title.x = element_blank())+
      scale_color_manual(values = c("#E41A1C","#377EB8"), name = "Treatment", labels = c("control", "solA")) +
      scale_shape_discrete(name = "Treatment", labels = c("control", "solA"))
  }
  return(g)
}
