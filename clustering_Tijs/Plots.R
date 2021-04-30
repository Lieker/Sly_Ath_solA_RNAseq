Avar <- function(genes, selection, samples){
  
  gene <- unlist(genes)
  samples <- as.data.frame(samples)
  selection <- as.matrix(selection)
  samples$gene <- row.names(samples)
  ConSam <- dcast(samples, condition ~ gene)
  ConSam <- as.data.frame(ConSam)
  row.names(ConSam) <- ConSam$condition
  ConSam$condition <- NULL
  ConSam <- as.matrix(ConSam)

  AV <- matrix(nrow = length(gene), ncol = nrow(ConSam))
  AV <- as.data.frame(AV, row.names = gene)
  colnames(AV) <- rownames(ConSam)
  
  ST <- matrix(nrow = length(gene), ncol = nrow(ConSam))
  ST <- as.data.frame(ST, row.names = gene)
  colnames(ST) <- rownames(ConSam)
  
  for(g in gene){
    avar <- list()
    stand <- list()
    for(i in row.names(ConSam)){
      l <- ConSam[i,]
      l <- as.vector(l)
      l <- l[!is.na(l)]
      avar[i] <- mean(selection[g,l])
      stand[i] <- sd(selection[g,l])
    }
    exp1 <- as.data.frame(avar)
    rownames(exp1) <- g
    AV[g,] <- exp1
    exp2 <- as.data.frame(stand)
    rownames(exp2) <- g
    ST[g,] <- exp2
  }
  AvCountsLong   <- melt(as.matrix(AV))
  StCountsLong   <- melt(as.matrix(ST))
  
outPut <- list()
outPut$AV <- AvCountsLong
outPut$ST <- StCountsLong

return(outPut)
}


exampleFeatureCounts <- function(){
  
  x <- data.frame("geneID"=c("gene1","gene2","gene3","gene4","gene5"),"Chr"=c("Chr1","Chr1","Chr1","Chr2","Chr2"),	
      "Start"=c("6860", "7678", "11192","0455", "1862"), "End"=c("7165","8121","9790", "1153","2567"),
      "Strand"=c("+","+","-","-","+"),	"Length"=c("269","444","536","699","706"),
    "sample1" = c("21","0","150", "12","0"),"sample2" = c("22","0","105","15","0"),"sample3" = c("31","1","165", "18","0"),"sample4" = c("2","0","83", "9","0"))
 
  return(x)
}

exampleProcessedCounts <- function(){
  
  x <- data.frame("geneTD"=c("gene1","gene2","gene3","gene4","gene5"),
                  "sample1" = c("21","0","150", "12","0"),"sample2" = c("22","0","105","15","0"),"sample3" = c("31","1","165", "18","0"),"sample4" = c("2","0","83", "9","0"))
  
  return(x)
}

exampleSample <- function(){
  
  x <- data.frame("sample"=c("sample1","sample2","sample3","sample4"),
                  "genotype" = c("wt","wt","mutant","mutant"),"condition" = c("control","treated","control","treated"))
  
  return(x)
}


