library(biomartr)
library(tidyverse)
library("biomaRt")
library(org.At.tair.db)
library(stringr)

get_Slycopersicum_annotations <- function(
) {
  descriptions <- read.csv("inputs/ITAG4.1_descriptions.csv", stringsAsFactors = FALSE)
  descriptions$gene <- substr(descriptions$gene, 1, 16)
  
  go_terms <- read.csv("inputs/ITAG4.1_goterms.csv", stringsAsFactors = FALSE)
  go_terms$gene <- substr(descriptions$gene, 1, 16)
  
  annotations <- left_join(descriptions, go_terms, by = "gene")
  return(annotations)
}
  



