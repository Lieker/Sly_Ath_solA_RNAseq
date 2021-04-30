introductie1<-function(){

return("
<br>
<h2> Introduction </h2><br>
<br>
RNA-Seq is a widely used method for studying the behavior of genes under different biological conditions.
To get a clear picture of what is happening in the dataset, it might be a valuable thing to do a 
cluster analysis. This will produce a nice picture of groups of genes that react(are regulated) in a similar way.<br>
In this application a number of different methods of clustering is combined to have a quick look wich fits your data best.
<br>
<br>
<h3>Input</h3><br>
For functioning the app needs two input files: 
<ol><li>A counts file as produced by <a href=\"http://subread.sourceforge.net/\">featureCounts</a> or a preprocessed file containing
only the read counts (tab delimitted, examples are depicted below.).</blockquote></li>
<li> A sample file (or design file) that should be a tab delimited file containing two or three columns with header.
First comlumn named <b><u>sample</u></b> containing the sample names <b><u>identical</u></b> in <b><u>name</u></b> and <b><u>order</u></b> as the sample names in the 
counts file followed by one or two columns containing conditions (genotypes, treatment, age, etc. depending on the experiment).</li></ol>
<br>
A featureCounts output file looks something like this:")
}

header1<-function(){
return("A preprocessed file should have the following lay-out:
")
}
header2<-function(){
return("The sample/experimental design file should look like this, including the samples and 1 or 2 columns of conditions (in this case genotype and treatment)
")
}


introductie2<-function(){
return("# Program:featureCounts v1.6.3; Command:\"featureCounts\" \"-O\" \"-t\" \"transcript\" \"-g\" \"gene_id\" \"-F\" \"gtf\" \"-a\" \"stringtie_transcriptome.gtf\" \"-o\"
 \"temp/results/counts.txt\"")
}

introductie3<-function(){
return("
Full example files to give the app a try, can be found <a href=\"https://github.com/KoesGroup/RNAseqCluster/tree/master/exampleFiles\">here</a>
<br>
<br>
<h3>Normalization</h3><br>
Normalization is done with <a href=\"http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html\">DESeq2</a>. 
All the different conditions (taken from the sample file) will be compared and genes that are differential in any of the 
possible combinations will be selected. This step does need a bit of time (30 sec. to 3 min.). After normalization a heatmap and a 
table containing the differientially expressed will be displayed, first 5 (head) or the full list (all).<br>
<br>
<br>
<h3>determination of the optimal cluster number</h3><br>
In the tab \"cluster number\" multiple algorithms are given to get an indication on how many clusters the data should be clustered in.
The outcome can be used for the clustering, but has to be seen as an indication.
In the end it is actualy up to the user to choose a number that fits the data best (or does best in grouping genes to your liking).<br>
<br>
<br>
<h3>The Clustering</h3><br>
Clustering can be done with multiple methods. It is difficult in advance to say what will work and what not, all methods have
there advantages and disadvantages, For more info on the different methods of clustering check the literature(google).<br>
As all these methods run pretty quick (depending on the total number of differiential genes) it is possible to run the
different methods with different settings and see what looks best. The clustering is visualized by a heatmap with a color bar indicating the cluster that the genes are clustered to,
and an expression plot of the cores/centers of clusters.<br>
<br>
<br>
<h3>Individual clusters</h3><br>
After the clustering you can select a cluster. After pressing submit a heatmap, expression plot and the list of genes in the cluster
will be displayed. If the clusters are not to your liking it is possible to go back to the clustering and rerun it with different settings.
<br>
<br>
<h3>Individual genes</h3><br>
To have a look at the expression if individual genes, write (or copy) the name of the gene you are interested and a barplot will be displayed.
If you want to have a look at more genes just add multiple in the box separated by a comma(without any spaces).<br>
<br>
<br>
<hr>
<br>
For questions or remarks, please feel free to send me an <a href=\"mailto:bliek@uva.nl SUBJECT=RNAseq_clustering\"> email.</a>
<br>
<br>
Tijs bliek,<br>
Universiteit van Amsterdam.<br>
<br>
Please check our <a href=\"https://www.scienceparkstudygroup.info\">webside: UvA study group</a>
<br>
<br>
<br>
<br>
")
}