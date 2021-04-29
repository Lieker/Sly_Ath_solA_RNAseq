# RNA-seq processed data

<!-- MarkdownTOC autolink="True" -->

- [File listing](#file-listing)
- [Raw RNA-seq fastq files](#raw-rna-seq-fastq-files)
- [Analytical pipeline](#analytical-pipeline)
- [Lexogen QuantSeq](#lexogen-quantseq)

<!-- /MarkdownTOC -->

This folder contains the processed data issued from the analysis of raw mRNA-seq from various genotypes affected in MYC1 expression. 

## File listing

- `config.yaml`: contains the settings used to run the pipeline. 
- `mapping_summary.tsv`: contains the statistics related to the RNA-seq alignment step (STAR software). Tab-separated file. 
- `multiqc_report.html`: contains the statistics related to read quality and adapter trimming. HTML file (open it in a web browser).
- `raw_counts.parsed.tsv`: the raw counts to be used for differential analysis. Tab-separated file. 
- `scaled_counts.tsv`: the scaled counts to be used for PCA, clustering, etc. __Do not use__ for differential analysis.  
- `samples.tsv`: the list of samples used. Feel free to add your experimental conditions (e.g. genotype) to perform the differential. analysis.  

## Raw RNA-seq fastq files

Ultimately, raw sequencing data have to be deposited at the [EBI European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home) before publication. 


## Analytical pipeline 
The raw fastq files were analysed using the following Snakemake-based pipeline, genomic references and parameters.  

- Snakemake mRNA-seq pipeline v0.3.5: [link](https://zenodo.org/record/4320989)
- Genomic references available [here](https://zenodo.org/record/4321000): 
    - Arabidopsis assembly: `GCA_000001735.2_TAIR10.1_genomic.fna`.
    - Arabidopsis genome annotation in GTF format is called `Araport11_GTF_genes_transposons.Mar92021.gtf`.
- Configuration and parameters used is visible in the [configuration file called `config.yaml`](./config.yaml)


## Lexogen QuantSeq
- Nature Methods paper: https://www.lexogen.com/wp-content/uploads/2015/04/nmeth.f.376.pdf
- Data analysis alternative STAR mapping parameters (more multimappers allowed): https://www.lexogen.com/quantseq-data-analysis-rev/ 
