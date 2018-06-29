# MOSCA

**Meta-Omics Software for Community Analysis**: a pipeline for performing Metagenomics and Metatranscriptomics integrated data analysis, in a local and fully automated workflow

## Features
* **Preprocessing** that starts with an initial quality check (FastQC), before Illumina artificial sequences removal and quality trimming (Trimmomatic, based on FastQC reports) and finally rRNA removal (SortMeRNA)
* **Assembly** that includes assembly with two optional tools (MetaSPAdes and Megahit) and quality control, with several metrics on the nature of contigs (MetaQUAST) and alignment of reads for estimating percentage of reads used in assembly (Bowtie2)
* **Annotation** that begins with gene calling (FragGeneScan), then annotation of identified ORFs (DIAMOND)
* **Analysis** where UniProt ID mapping is performed (if UniProt database was chosen) for retrieving taxonomic and functional information as well as, if MT data was included, differential gene expression and multisample comparison (HTSeq-count and DeSEQ2)

## Setting up MOSCA

MOSCA files can be retrieved with a git command.

```
git clone https://github.com/iquasere/MOSCA.git
```

MOSCA already brings a bash script that will install all of its pre-dependencies and databases.

```
cd MOSCA
bash install.bash
```

## Base arguments for running MOSCA

MOSCA was designed to run with as few arguments as possible. Only the input files and output directory have to be specified. MOSCA considers the input files in the following format, separated by spaces:

```
path/to/mg_file1,path/to/mg_file2:path/to/mt_file1,path/to/mt_file2
```

MOSCA is prepared to handle experiments input if some of these files is not available (because the experiment was single-end or metatranscriptomic data was not analyzed). How to format such cases is explained in the next examples. Note that this are the examples for running MOSCA from its directory ('MOSCA'). To run them from other directory, the path to the main script ('main.py') must be specified.

* One single-end file

```
python main.py --files path/to/file --output-dir output_directory
```

* Two single-end files (two different experiments)

```
python main.py --files path/to/file1 path/to/file2 --output-dir output_directory
```

* Two paired-end files

```
python main.py --files path/to/file1,path/to/file2 --output-dir output_directory
```

* Four paired-end files with MT data (two different experiments)

```
python main.py --files path/to/mg_file1_of_experiment1,path/to/mg_file2_of_experiment1:path/to/mt_file1_of_experiment1,path/to/mt_file2_of_experiment1 path/to/mg_file1_of_experiment2,path/to/mg_file2_of_experiment2:path/to/mt_file1_of_experiment2,path/to/mt_file2_of_experiment2 --output-dir output_directory
```
