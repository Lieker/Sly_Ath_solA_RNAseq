setwd("C:/Users/levla/github/Ath_RNAseq/Sly1")
source("scripts/volcanoplot.R")
source("scripts/plsda.R")
source("scripts/permanova.R")
source("scripts/create_pca_plot.R")
source("scripts/GO_annotation_list.R")
source("scripts/get_annotated_DEGs.R")

make_volcanoplot(counts_csv_file = "inputs/raw_counts.csv",
                 xp_design_csv_file = "inputs/xp_design.csv",
                 plantpart = "root",
                 ref_treatment = "no_solA",
                 treatment2 = "yes_solA",
                 log2FC_threshold = 0,
                 FCcutoff_volcano = 1,
                 padj_threshold = 0.05)
do_plsda(xp_design_csv_file = "inputs/xp_design.csv",
         counts_csv_file = "inputs/raw_counts.csv",
         organ = "root",
         nr = 9)
do_permanova(xp_design_csv_file = "inputs/xp_design.csv",
             counts_csv_file = "inputs/raw_counts.csv",
             plantpart = "root")
plot_pca(count_csv_file = "inputs/raw_counts.csv", 
         xp_design_csv_file = "inputs/xp_design.csv", 
         plantpart = "root", 
         pc_x_axis = "PC1", 
         pc_y_axis = "PC2", 
         pca_colour = "treatment")
GO_annotation_list(counts_csv_file = "inputs/raw_counts.csv",
                   xp_design_csv_file = "inputs/xp_design.csv",
                   timepoint = 2,
                   ref_treatment = "ethanol",
                   treatment2 = "millimolar_solanoeclepinA",
                   log2FC_threshold = 0,
                   padj_threshold = 0.05,
                   organism = "Arabidopsis thaliana",
                   attr = c("tair_symbol", 
                            "uniprotswissprot",
                            "entrezgene_id",
                            "description",
                            "external_gene_name",
                            "external_gene_source"),
                   name = "annotated_DEGslist.csv",
                   ont = "BP", #can be either BP, CC or MF
                   namego = "outputs/GOannotationslist.csv")
get_annotated_DEGs <- function(counts_csv_file = "inputs/raw_counts.csv",
                               xp_design_csv_file = "inputs/xp_design.csv",
                               plantpart = "root",
                               ref_treatment = "no_solA",
                               treatment2 = "yes_solA",
                               log2FC_threshold = 0,
                               padj_threshold = 0.05,
                               organism = "Solanum lycopersicum",
                               attr = c("description",
                                        "athaliana_eg_homolog_ensembl_gene",
                                        "athaliana_eg_homolog_associated_gene_name",
                                        "external_gene_name"),
                               name = "outputs/annotated_DEGslist.csv")