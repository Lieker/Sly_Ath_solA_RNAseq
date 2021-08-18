
produce_GOdf <- function(r, d,
                         GOid = c("GO:0006952","GO:0045087")){
  GOdf <- res_annotated[res_annotated$go_id %in% GOid,]
  uniquepaste <- function(x) { paste(unique(x), sep = ',', collapse = ",") }
  GOdf <- aggregate(GOdf, by = list(GOdf$transcript), FUN = uniquepaste)
  t <- colData(d)
  for(i in GOdf$transcript){
    d <- plotCounts(gene = i, 
                    dds = dds_all, 
                    intgroup = "treatment", 
                    returnData = TRUE)
    t$c <- d$count
    colnames(t)[which(names(t) == "c")] <- paste0("counts_",i)
  }
  t <- melt(as.data.frame(t), id.vars = c("sample","treatment","solA","N","P","sizeFactor"))
  return(t)
}
