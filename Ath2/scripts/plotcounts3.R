plotcounts3 <- function(GOdf, d = dds, GOlist) {
  GOdf <- GOdf[GOdf$go_id %in% GOlist,]
  uniquepaste <- function(x) { paste(unique(x), sep = ',', collapse = ",") }
  GOdf <- aggregate(GOdf, by = list(GOdf$transcript), FUN = uniquepaste)
  
  g <- ggplot()+
    scale_color_brewer(palette="Set1") +
    theme_classic() +
    scale_y_log10(breaks=c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000))
  
  for(i in GOdf$transcript){
    d <- plotCounts(d, gene = i,
                    intgroup = c("treatment","solA","N","P"), 
                    returnData = TRUE) %>% rownames_to_column("sample")
    
    colnames(d)[2] <- "counts"
    
    d$solA <- gsub("yes", "SolA", d$solA)
    d$solA <- gsub("no", "cntrl", d$solA)
    d <- data.frame(lapply(d, function(x) {gsub("yes", "+", x)}))
    d <- data.frame(lapply(d, function(x) {gsub("no", "-", x)}))
    
    d$group <- paste("N",d$N,"P",data$P)
    d <- data.frame(lapply(d, function(x) {gsub("N ", "N", x)}))
    d <- data.frame(lapply(d, function(x) {gsub("P ", "P", x)}))
    
    d$solA <- factor(d$solA, levels = c("cntrl","SolA"))
    d$group <- as.factor(d$group)
    d$treatment <- as.factor(d$treatment)
    d$counts <- as.numeric(as.character(d$counts))
    data.summary <- d %>% group_by(treatment,group,solA) %>% summarise(sd = sd(counts, 
                                                                               na.rm = TRUE),
                                                                       meancounts = mean(counts))
    g <- g +
      geom_line(data = data.summary,
                aes(x = group,
                    y = meancounts,
                    colour = solA,
                    group = solA),
                size = 1.2) +
      geom_point(data = data.summary,
                 aes(x = group,
                     y = meancounts),
                 size = 3) 
  }
  return(g)
}
