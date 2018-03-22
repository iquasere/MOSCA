# MOSCA

**Meta-Omics Software for Community Analysis**: a pipeline for performing Metagenomics and Metatranscriptomics integrated data analysis, in a local and fully automated workflow

## Features
* **Preprocessing** that starts with an initial quality check (FastQC), before Illumina artificial sequences removal and quality trimming (Trimmomatic, based on FastQC reports) and finally rRNA removal (SortMeRNA)
* **Assembly** that includes assembly with two optional tools (MetaSPAdes and Megahit) and quality control, with several metrics on the nature of contigs (MetaQUAST) and alignment of reads for estimating percentage of reads used in assembly (Bowtie2)
* **Annotation** that begins with gene calling (FragGeneScan), then annotation of identified ORFs (DIAMOND)
* **Analysis** where UniProt ID mapping is performed (if UniProt database was chosen) for retrieving taxonomic and functional information as well as, if MT data was included, differential gene expression and multisample comparison (HTSeq-count and DeSEQ2)
