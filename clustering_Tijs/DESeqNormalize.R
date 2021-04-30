normalize <- function(countdata, samplefile, f, p){
  col.names = colnames(samplefile)                                                # extract all column names
  columns.to.discard = c("sample","fq1","fq2","fq")                                    # column we don't need to extract all conditions
  colsForConditions = col.names[! col.names %in% columns.to.discard]              # only keeping the column of interest
    # one condition
  if (length(colsForConditions) == 1){
    condition <- factor(samplefile[,colsForConditions])
    # two conditions
  } else if (length(colsForConditions == 2)){
    # two conditions --> make a third column that is the product of the two
    samplefile$conditions = paste(samplefile[,colsForConditions[1]],samplefile[,colsForConditions[2]],sep = ".")
    condition <- factor(x = samplefile[,"conditions"],levels = unique(samplefile[,"conditions"]))
  } else if (length(conditions > 2)){
    print("too many conditions to compare, skipping")
  }
  helperFile <- as.data.frame(condition)
  rownames(helperFile) <- samplefile$sample
  #colnames(helperFile) <- NULL
  # Analysis with DESeq2 ----------------------------------------------------
  # Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  coldata <- data.frame(row.names=colnames(countdata), condition)
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
  # Run the DESeq pipeline
  dds <- DESeq(dds)
  # create dataframe containing normalized counts, to which the differential expression values will be added
  resdata <- as.data.frame(counts(dds, normalized=TRUE))
  #p <- 0.01
  #f <- 2.0
  # iterate through the list f conditions to create differential expression (DE) values for all possible pairs
  total_selec <- list()
  x <- 1
  for(i in levels(condition)){
    x <- x + 1
    if (x <= length(levels(condition))){
      for(j in x:length(levels(condition))){
        res <- results(dds, contrast=c("condition", i,  levels(condition)[j]))                      # i = first in pair, levels(condition)[j] is the second in pair.
        d <- paste(i, levels(condition)[j], sep="&")                                                # paste the two conditions in one character, to be used as the pair name
        res$genenames <- rownames(res)
        resul <- as.data.frame(res)
        goede <- resul %>% filter(padj<p & (log2FoldChange>f | log2FoldChange<(-1*f)) )
        selec <- as.list(goede$genenames)
        total_selec <- append(total_selec, selec)
        print(length(selec))# get number of DE values with P-value < 0.05
      }
    }
  }
  total_selec <- c(unique(total_selec))
  total_selec <- t(as.data.frame(total_selec))
  selection <- resdata[total_selec[,1],]
  
  samples <- as.data.frame(helperFile)
  samples$gene <- row.names(samples)
  ConSam <- dcast(samples, condition ~ gene)
  ConSam <- as.data.frame(ConSam)
  row.names(ConSam) <- ConSam$condition
  ConSam$condition <- NULL
  ConSam <- as.matrix(ConSam)
  scaledata <- t(scale(t(selection))) # Centers and scales data.
  scaledata <- scaledata[complete.cases(scaledata),]
  scaledata1 <- as.data.frame(scaledata)
  mydat <- data.frame(row.names = rownames(scaledata))
  
  
  for(i in row.names(ConSam)){
    l <- ConSam[i,]
    l <- as.vector(l)
    l <- l[!is.na(l)]
    mydat[i] <- (rowMeans(as.data.frame(scaledata1[unlist(c(l))])))
  }
  mydat <- as.matrix(mydat)  
  

#  samples$gene <- row.names(samples)
#  ConSam <- dcast(samples, condition ~ gene)
#  ConSam <- as.data.frame(list(split(apply(ConSam, 1, function(x) unname(x[!is.na(x)])), 1:nrow(ConSam))),, row.names = 1)
  
#  ConSam <- as.matrix(ConSam)
#  scaledata <- t(scale(t(selection))) # Centers and scales data.
#  scaledata <- scaledata[complete.cases(scaledata),]
#  scaledata1 <- as.data.frame(scaledata)
#  mydat <- data.frame(row.names = rownames(scaledata))
#  
#  for(i in rownames(ConSam)){
#    l <- ConSam[i,]
#    l <- as.vector(l)
#    names(l) <- NULL
#    mydat[i] <- (rowMeans(as.data.frame(scaledata1[unlist(c(l))])))
#  }
#  mydat <- as.matrix(mydat)
  
  
  
  
  normData <- list()
  normData$samples <- helperFile
  normData$selection <- selection
  normData$scaledata <- scaledata
  normData$avarages <- mydat
  return(normData)
}
