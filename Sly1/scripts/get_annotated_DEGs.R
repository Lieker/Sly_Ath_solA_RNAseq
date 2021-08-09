source("Sly1/scripts/get_filtered_list_of_DEGs.R")
source("Sly1/scripts/get_Slycopersicum_annotations.R")

get_annotated_DEGs <- function(counts_csv_file = "Sly1/input/counts.csv",
                               xp_design_csv_file = "Sly1/input/xp_design.csv",
                               plantpart = "root",
                               ref_treatment = "no_solA",
                               treatment2 = "yes_solA",
                               log2FC_threshold = 0,
                               padj_threshold = 0.05,
                               name = "Sly1/output/annotated_DEGslist.csv",
                               attr = c("description",
                                        "athaliana_eg_homolog_ensembl_gene",
                                        "athaliana_eg_homolog_associated_gene_name",
                                        "external_gene_name",
                                        "athaliana_eg_homolog_perc_id"),
                               method
                               ) {
  res <- get_list_of_DEGs(counts_csv_file, 
                         xp_design_csv_file, 
                         plantpart, 
                         ref_treatment, 
                         treatment2, 
                         log2FC_threshold, 
                         padj_threshold,
                         method) %>% rownames_to_column("original_genename")
  genes <- res$original_genename
  genes <- substr(genes, start = 6, stop = 21)
  res$gene <- genes
     
  annotations <- get_Slycopersicum_annotations(attr = attr)
  res_annotated <- left_join(res, annotations, by = "gene")
  return(res_annotated)
  
  uniquepaste <- function(x) { paste(unique(x), sep = ',', collapse = ",") }
  res_annotated <- aggregate(res_annotated, by = list(res_annotated$gene), FUN = uniquepaste)
  write_delim(res_annotated, file = name, delim = ",", col_names = TRUE)
  
  return(res_annotated)
}
