![ScreenShot](mosca_logo.png)

Logo by Sérgio A. Silva

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

## MOSCA is finally available as a Docker image!

To use MOSCA's Docker version, Docker must first be installed.

```
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

After having docker installed in the system, pull MOSCA's image.

```
docker pull iquasere/mosca:latest
```

At this point, MOSCA allows for defining custom databases for annotation, but not for adapters nor rRNA identification (instead, the image already brings the databases from Trimmomatic and SortMeRNA, respectively).
The database(s) for annotation must be present in FASTA or DMND (diamond binary) format in a specific directory that must be referenced in the mosca command.

```
docker run -it -v /path/to/folder_of_databases:/MOSCA/Databases/annotation_databases iquasere/mosca [arguments]
```

"/path/to/folder_of_databases" is the directory where the databases are stored. The rest of the command is to be inputed exactly as presented here, except for the [arguments], which are to be inputed just like if not using docker.

## Base arguments for running MOSCA

MOSCA was designed to run with as few arguments as possible. Only the input files and output directory have to be specified. MOSCA considers the input files in the following format, separated by spaces:

```
path/to/mg_file1,path/to/mg_file2:path/to/mt_file1,path/to/mt_file2
```

MOSCA is prepared to handle experiments input if some of these files is not available (because the experiment was single-end or metatranscriptomic data was not analyzed). How to format such cases is explained in the next examples. Note that this are the examples for running MOSCA from its directory ('MOSCA'). To run them from other directory, the path to the main script ('main.py') must be specified.

* One single-end file

```
python MOSCA/mosca.py --files path/to/file --output-dir output_directory
```

* Two single-end files (two different experiments)

```
python MOSCA/mosca.py --files path/to/file1 path/to/file2 --output-dir output_directory
```

* Two paired-end files

```
python MOSCA/mosca.py --files path/to/file1,path/to/file2 --output-dir output_directory
```

* Four paired-end files with MT data (two different experiments)

```
python MOSCA/mosca.py --files path/to/mg_file1_of_experiment1,path/to/mg_file2_of_experiment1:path/to/mt_file1_of_experiment1,path/to/mt_file2_of_experiment1 path/to/mg_file1_of_experiment2,path/to/mg_file2_of_experiment2:path/to/mt_file1_of_experiment2,path/to/mt_file2_of_experiment2 --output-dir output_directory
```


## Additional parameters

MOSCA includes the option to chose between MetaSPAdes (default) and Megahit as the tool for assembly. Also, while many functions are developed only for UniProt, another database can be chosen. Since not all steps are always wanted, MOSCA allows to chose whether or not to perform every step of its analysis.

```
usage: mosca.py [-h] [-f [FILES [FILES ...]]] [-data {paired,single}]
                [-a Assembler] [-db Database] [-o Directory] [-nopp] [-noas]
                [-noan] [-nobin] [-ol {min,med,max}] [-mp]
                [-c [CONDITIONS [CONDITIONS ...]]]

Multi Omics Software for Community Analysis

optional arguments:
  -h, --help            show this help message and exit
  -f [FILES [FILES ...]], --files [FILES [FILES ...]]
                        Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2
                        mg2R1,...)
  -data {paired,single}, --type-of-data {paired,single}
                        Type of data (paired/single)-ended
  -a Assembler, --assembler Assembler
                        Tool for assembling the reads
  -db Database, --annotation-database Database
                        Database for annotation (.fasta or .dmnd)
  -o Directory, --output Directory
                        Directory for storing the results
  -nopp, --no-preprocessing
                        Don't perform preprocessing
  -noas, --no-assembly  Don't perform assembly
  -noan, --no-annotation
                        Don't perform annotation
  -nobin, --no-binning  Don't perform binning
  -ol {min,med,max}, --output-level {min,med,max}
                        Level of file output from MOSCA, min outputs only the
                        analysis results, med removes intermediate files, max
                        outputs all intermediate and final data
  -mp, --metaproteomic  If data is metagenomic and metaproteomic, if not
                        specified will be assumed to be metagenomic and
                        metatranscriptomic
  -c [CONDITIONS [CONDITIONS ...]], --conditions [CONDITIONS [CONDITIONS ...]]
                        Different conditions for metatranscriptomic analysis

A tool for performing metagenomics, metatranscriptomics and metaproteomics
analysis.
```
