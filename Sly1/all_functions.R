setwd("C:/Users/levla/github/Ath_RNAseq/Sly1")
source("Sly1/scripts/volcanoplot.R")
source("Sly1/scripts/plsda.R")
source("Sly1/scripts/permanova.R")
source("Sly1/scripts/create_pca_plot.R")
source("Sly1/scripts/GO_annotation_list.R")
source("Sly1/scripts/get_annotated_DEGs.R")

make_volcanoplot(counts_csv_file = "inputs/raw_counts.csv",
                 xp_design_csv_file = "inputs/xp_design.csv",
                 timepoint = 2,
                 ref_treatment = "ethanol",
                 treatment2 = "millimolar_solanoeclepinA",
                 log2FC_threshold = 0,
                 FCcutoff_volcano = 0,
                 padj_threshold = 0)
do_plsda(xp_design_csv_file = "inputs/xp_design.csv",
         counts_csv_file = "inputs/raw_counts.csv",
         timepoint = 2,
         nr = 999)
do_permanova(xp_design_csv_file = "inputs/xp_design.csv",
             counts_csv_file = "inputs/raw_counts.csv",
             timepoint = 2)
plot_pca(count_csv_file = "inputs/raw_counts.csv", 
         xp_design_csv_file = "inputs/xp_design.csv", 
         timepoint = 2, 
         pc_x_axis = "PC1", 
         pc_y_axis = "PC2", 
         pca_colour = "sequencing_batch")
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
get_annotated_DEGs(counts_csv_file = "inputs/raw_counts.csv",
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
                   name = "outputs/annotated_DEGslist.csv")


#try to make 3 lists for each timeponit/comparison. em2 is taking into account both rna_isolation_batch
#and purification_procedure, em2_ is only corrected for rna_isolation_Batch, em2__ for purification_
#procedure. So change DEG_function function every time and repeat this 2x.

em2__ <- get_annotated_DEGs()
em6__ <- get_annotated_DEGs(timepoint = 6)
em24__ <- get_annotated_DEGs(timepoint = 24)
en2__ <- get_annotated_DEGs(treatment2 = "nanomolar_solanoeclepinA")
en6__ <- get_annotated_DEGs(timepoint = 6, treatment2 = "nanomolar_solanoeclepinA")
en24__ <- get_annotated_DEGs(timepoint = 24, treatment2 = "nanomolar_solanoeclepinA")
nm2__ <- get_annotated_DEGs(ref_treatment = "nanomolar_solanoeclepinA", treatment2 = "millimolar_solanoeclepinA")
nm6__ <- get_annotated_DEGs(timepoint = 6, ref_treatment = "nanomolar_solanoeclepinA", treatment2 = "millimolar_solanoeclepinA")
nm24__ <- get_annotated_DEGs(timepoint = 24, ref_treatment = "nanomolar_solanoeclepinA", treatment2 = "millimolar_solanoeclepinA")

em2$batchcorr <- "r+p"
nm6$batchcorr <- "r+p"
em2_$batchcorr <- "r"
em24_$batchcorr <- "r"
em6_$batchcorr <- "r"
en2_$batchcorr <- "r"
en6_$batchcorr <- "r"
nm2_$batchcorr <- "r"
nm24_$batchcorr <- "r"
nm6_$batchcorr <- "r"
em2__$batchcorr <- "p"
em24__$batchcorr <- "p"
en2__$batchcorr <- "p"
en24__$batchcorr <- "p"
en6__$batchcorr <- "p"
nm2__$batchcorr <- "p"
nm24__$batchcorr <- "p"

em2$refgroup <- "em2"
nm6$refgroup <- "nm6"
em2_$refgroup <- "em2"
em24_$refgroup <- "em24"
em6_$refgroup <- "em6"
en2_$refgroup <- "en2"
en6_$refgroup <- "en6"
nm2_$refgroup <- "nm2"
nm24_$refgroup <- "nm24"
nm6_$refgroup <- "nm6"
em2__$refgroup <- "em2"
em24__$refgroup <- "em24"
en2__$refgroup <- "en2"
en24__$refgroup <- "en24"
en6__$refgroup <- "en6"
nm2__$refgroup <- "nm2"
nm24__$refgroup <- "nm24"

alldegs <- rbind(em2, em2_, em2__, em24_, em24__, em6_, en2_, en2__, en24__, en6_, en6__, nm2_, nm2__, nm24_, nm24__, nm6, nm6_)

all_arabidopsis_genes_annotated <- get_Athaliana_annotations(attr = c("affy_ath1_121501", "entrezgene_id"))

degs_n_distinct <- degs_n %>% distinct(genes, baseMean, batchcorr, refgroup, .keep_all = TRUE)

write.csv(degs_n_distinct, "outputs/alldegs.csv", row.names = F)


