![ScreenShot](mosca_logo.png)

Logo by [Sérgio A. Silva](https://www.ceb.uminho.pt/People/Details/64888072-5cde-42da-b7e5-691d380cefb2)

# Meta-Omics Software for Community Analysis

**MOSCA** (portuguese for fly) is a pipeline designed for performing Metagenomics and Metatranscriptomics integrated data analysis, in a mostly local and fully automated workflow.

## Features
* **Preprocessing** where low quality regions of data are trimmed and reads less interest are removed. FastQC's reports are used to automatically set the parameters for the other tools. It includes:
    * initial quality check with **FastQC**
    * Illumina artificial sequences removal with **Trimmomatic**: based on **FastQC** reports, MOSCA will find the adapters file most approprita to the data
    * rRNA removal with **SortMeRNA**: uses Pfam and SILVA databases as reference
    * quality trimming with **Trimmomatic**: 
        * another **FastQC** report will be generated after rRNA removal, and will be used to set the parameters for **Trimmomatic**'s hard trimmers (CROP and HEADCROP). This will ensure that the data will be reported as excellent by FastQC
        * reads with less than 20 average quality or 100 nuleotides of length will also be removed
    * final quality check with **FastQC** 
* **Assembly** where metagenomics (MG) trimmed reads will be assembled to partially reconstruct the original genomes in the samples. It includes:
    * assembly with two possible assemblers - **MetaSPAdes** and **Megahit** - which will be used in a multi-kmer approach
    * control over the quality of the contigs, with **MetaQUAST** reporting on several classical metrics (such as N50 and L50) and alignment of reads for estimating percentage of reads used in assembly, with **Bowtie2**
* **Annotation** where proteins present in the contigs will be identified. It includes:
    * gene calling with **FragGeneScan**
    * annotation of identified ORFs with **DIAMOND**, using the **UniProt database** as reference - MOSCA only reports on the first annotation
    * retrieval of biological information with **UniProt ID mapping** API
    * functional annotation with **Reverse PSI-BLAST** (RPSBLAST), using the **COG database** as reference
        * MOSCA automatically **generates new databases by the number of threads specified**, thus allowing for multithread annotation with **RPSBLAST**
    * the quantification of each protein in MG data, by alignment of MG reads to the contigs using **Bowtie2** and quantification of reads to protein using **HTSeq-count**
* **MetaTranscriptomics (MT) analysis** where the expression of each identified protein is quantified. It includes:
    * alignment of MT reads to the MG contigs with **Bowtie2**, and quantification of reads to protein using **HTSeq-count**
    * differential gene expression and multisample comparison using **DeSEQ2**
* **Normalization** of protein quantification for the final report using **edgeR**

## Setting up MOSCA

MOSCA files can be retrieved with a git command.

```
git clone https://github.com/iquasere/MOSCA.git
```

MOSCA already brings a bash script that will install all of its pre-dependencies and databases. However, it does require Conda previously installed. Instructions on installing Anaconda on an Ubuntu 18.04 may be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart).

```
bash MOSCA/install.bash
```
<!---
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
-->
## Base arguments for running MOSCA

MOSCA was designed to run with as few arguments as possible. Only the input files and output directory have to be specified. MOSCA considers the input files in the following format, separated by spaces:

```
path/to/mg_file1,path/to/mg_file2:path/to/mt_file1,path/to/mt_file2
```

MOSCA is prepared to handle experiments input if some of these files is not available (because the experiment was single-end or metatranscriptomic data was not analyzed). How to format such cases is explained in the next examples. Note that this are the examples for running MOSCA from its directory ('MOSCA'). To run them from other directory, the path to the main script ('main.py') must be specified.

* One single-end file

```
python MOSCA/scripts/mosca.py --files path/to/file --output-dir output_directory
```

* Two single-end files (two different experiments)

```
python MOSCA/scripts/mosca.py --files path/to/file1 path/to/file2 --output-dir output_directory
```

* Two paired-end files

```
python MOSCA/scripts/mosca.py --files path/to/file1,path/to/file2 --output-dir output_directory
```

* Four paired-end files with MT data (two different experiments)

```
python MOSCA/scripts/mosca.py --files path/to/mg_file1_of_experiment1,path/to/mg_file2_of_experiment1:path/to/mt_file1_of_experiment1,path/to/mt_file2_of_experiment1 path/to/mg_file1_of_experiment2,path/to/mg_file2_of_experiment2:path/to/mt_file1_of_experiment2,path/to/mt_file2_of_experiment2 --output-dir output_directory
```


## Additional parameters

MOSCA includes the option to chose between MetaSPAdes (default) and Megahit as the tool for assembly. Also, while many functions are developed only for UniProt, another database can be chosen. Since not all steps are always wanted, MOSCA allows to chose whether or not to perform every step of its analysis.

```
usage: mosca.py [-h] [-f [Input files [Input files ...]]]
                [-data {paired,single}] [-a Assembler] [-db Database]
                [-o Directory] [-nopp] [-noas] [-noan] [-nobin]
                [-ol {min,med,max}] [-mp] [-c [CONDITIONS [CONDITIONS ...]]]
                [-t Threads] [-m Memory]

Multi Omics Software for Community Analysis

optional arguments:
  -h, --help            show this help message and exit
  -f [Input files [Input files ...]], --files [Input files [Input files ...]]
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
                        Different conditions for metatranscriptomic analysis,
                        separated by comma (,)
  -t Threads, --threads Threads
                        Number of threads available for MOSCA
  -m Memory, --memory Memory
                        Maximum memory (byte) available for MOSCA. Applied
                        only in the assembly

A tool for performing metagenomics, metatranscriptomics and metaproteomics
analysis.
```
