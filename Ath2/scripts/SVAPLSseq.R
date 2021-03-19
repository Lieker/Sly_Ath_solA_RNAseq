library(SVAPLSseq)
library(SummarizedExperiment)
library(edgeR)

source("scripts/filter_raw_counts.R")

c <- filter_raw_counts(trtm = c("a","b"))
rownames(c) <- NULL
class(c) <- "numeric"
c.se <- SummarizedExperiment(assays = SimpleList(counts = c))
c.se <- DGEList(counts = c)
xd <- as.factor(xp_design[1:10,2])
sv <- svplsSurr(dat = c.se, group = xd, max.surrs = 3, surr.select = "automatic", controls = NULL)