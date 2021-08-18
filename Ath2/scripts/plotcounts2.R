plotcounts2 <- function(gene, d = dds) {
  data <- plotCounts(d, gene,
                     intgroup = c("treatment","solA","N","P"), 
                     returnData = TRUE) %>% rownames_to_column("sample")
  
  
  colnames(data)[2] <- "counts"
  
  data$solA <- gsub("yes", "SolA", data$solA)
  data$solA <- gsub("no", "cntrl", data$solA)
  data <- data.frame(lapply(data, function(x) {gsub("yes", "+", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("no", "-", x)}))
  
  data$group <- paste("N",data$N,"P",data$P)
  data <- data.frame(lapply(data, function(x) {gsub("N ", "N", x)}))
  data <- data.frame(lapply(data, function(x) {gsub("P ", "P", x)}))
  
  data$solA <- factor(data$solA, levels = c("cntrl","SolA"))
  data$group <- as.factor(data$group)
  data$treatment <- as.factor(data$treatment)
  data$counts <- as.numeric(as.character(data$counts))
  data.summary <- data %>% group_by(treatment,group,solA) %>% summarise(sd = sd(counts, 
                                                                     na.rm = TRUE),
                                                             meancounts = mean(counts))
  g <- ggplot() + 
    geom_jitter(data = data, 
                aes(x = group,
                    y = counts,
                    colour = solA),
                position= position_jitter(0.2),
                alpha = 0.8) +
    geom_line(data = data.summary,
              aes(x = group,
                  y = meancounts,
                  colour = solA,
                  group = solA),
              size = 1.2)+
    scale_color_brewer(palette="Set1") +
    geom_point(data = data.summary,
               aes(x = group,
                   y = meancounts),
               size = 3) +
    theme_classic() +
    ggtitle(paste0(gene)) +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000))
  
  return(g)
}