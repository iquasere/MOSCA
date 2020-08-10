![ScreenShot](mosca_logo.png)

Logo by [Sérgio A. Silva](https://www.ceb.uminho.pt/People/Details/64888072-5cde-42da-b7e5-691d380cefb2)

# Meta-Omics Software for Community Analysis (MOSCA)

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
* **Binning** where the contigs are clustered into taxonomic units, to validate (or not) the annotation, and possibily help reconstructing genomes from the samples
    * **MaxBin2** bins the contigs by tetranucleotide composition, relative abundance, and marker genes analysis
    * the final bins are reported for their completeness - how many of the marker genes are present in each bin?
* **MetaTranscriptomics (MT) analysis** where the expression of each identified protein is quantified. It includes:
    * alignment of MT reads to the MG contigs with **Bowtie2**, and quantification of reads to protein using **HTSeq-count**
    * differential gene expression and multisample comparison using **DESeq2**
* **Normalization** of protein quantification for the final report using **edgeR**
* **Pathway representation** in **KEGG maps**, representing both the metabolic networks of most abundant taxa and expression levels of metabolic functions

## Setting up MOSCA

MOSCA files can be retrieved from git.
```
git clone https://github.com/iquasere/MOSCA.git
```

MOSCA requires Conda previously installed. Instructions on installing Anaconda on an Ubuntu 18.04 may be found [here](https://www.digitalocean.com/community/tutorials/how-to-install-anaconda-on-ubuntu-18-04-quickstart). Alternatively, you can directly download the installation for Ubuntu with which MOSCA was tested - problems have been found using more recent distributions.
```
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
bash Anaconda3-2019.03-Linux-x86_64.sh
```

After having Conda available, the bash script that will install all of its pre-dependencies and databases can be run as follow: 
```
bash MOSCA/install.bash
```

## Base arguments for running MOSCA

MOSCA was designed to run with as few arguments as possible. Only the input files and output directory have to be specified. MOSCA considers the input files in the following format, separated by spaces:
```
path/to/mg_file1,path/to/mg_file2:path/to/mt_file1,path/to/mt_file2
```

MOSCA is prepared to handle experiments input if some of these files is not available (because the experiment was single-end or metatranscriptomic data was not analyzed). How to format such cases is explained in the next examples. Note that this are the examples for running MOSCA from its directory ('MOSCA'). To run them from other directory, the path to the main script ('main.py') must be specified.

* One single-end file

```
python MOSCA/scripts/mosca.py --files path/to/file --output output_directory
```

* Two single-end files (two different experiments)

```
python MOSCA/scripts/mosca.py --files path/to/file1 path/to/file2 --output output_directory
```

* Two paired-end files

```
python MOSCA/scripts/mosca.py --files path/to/file1,path/to/file2 --output output_directory
```

* Four paired-end files with MT data (two different experiments)

```
python MOSCA/scripts/mosca.py --files path/to/mg_file1_of_experiment1,path/to/mg_file2_of_experiment1:path/to/mt_file1_of_experiment1,path/to/mt_file2_of_experiment1 path/to/mg_file1_of_experiment2,path/to/mg_file2_of_experiment2:path/to/mt_file1_of_experiment2,path/to/mt_file2_of_experiment2 --output output_directory
```


## Additional parameters

MOSCA includes the option to chose between MetaSPAdes (default) and Megahit as the tool for assembly. Also, while many functions are developed only for UniProt, another database can be chosen. Since not all steps are always wanted, MOSCA allows to chose whether or not to perform every step of its analysis.

```
usage: mosca.py [-h] [-f [Input files [Input files ...]]]
                [-st {paired,single}] [-a Assembler] [-db Database]
                [-o Directory] [-nopp] [-noas] [-noan] [-nobin]
                [-ol {min,med,max}]
                [-tod {metatranscriptomics,metaproteomics}]
                [-c [CONDITIONS [CONDITIONS ...]]] [-t Threads] [-m Memory]
                [-mark Marker gene set]

Meta-Omics Software for Community Analysis

optional arguments:
  -h, --help            show this help message and exit
  -f [Input files [Input files ...]], --files [Input files [Input files ...]]
                        Input files for analysis (mg1R1,mg1R2:mt1R1,mt1R2
                        mg2R1,...)
  -st {paired,single}, --sequencing-technology {paired,single}
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
  -tod {metatranscriptomics,metaproteomics}, --type-of-data {metatranscriptomics,metaproteomics}
                        If data is metagenomics integrated with
                        metatranscriptomics or metaproteomics, if not
                        specified will be assumed to be metagenomics and
                        metatranscriptomics. This option can be ignored with
                        dealing only with metagenomics data
  -c [CONDITIONS [CONDITIONS ...]], --conditions [CONDITIONS [CONDITIONS ...]]
                        Different conditions for
                        metatranscriptomics/metaproteomics analysis, separated
                        by comma (,)
  -t Threads, --threads Threads
                        Number of threads available for MOSCA
  -m Memory, --memory Memory
                        Maximum memory (byte) available for MOSCA. Applied
                        only in the assembly
  -mark Marker gene set, --marker-gene-set Marker gene set
                        Marker gene set to use for binning with MaxBin2. 107
                        if archaea are not to be considered, 40 if data is
                        diverse.

A tool for performing metagenomics, metatranscriptomics and metaproteomics
analysis.
```

## MOSCA is now also available through Bioconda!

To use MOSCA through Bioconda, an environment must be created containing all the dependencies. 
Create the environment and install MOSCA with ```conda create -n mosca -c bioconda -c conda-forge mosca```
Activate the environment with ```conda activate mosca```
Run MOSCA with ```mosca.py [arguments]```

## MOSCA is also available as a Docker image!

To use MOSCA's Docker version, Docker must first be installed.

```
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

After having docker installed in the system, pull MOSCA's image.

```
docker pull iquasere/mosca:latest
```

MOSCA's docker version allows for defining custom databases for annotation, but not for adapters nor rRNA identification (instead, the image already brings the databases from Trimmomatic and SortMeRNA, respectively).
The database(s) for annotation must be present in FASTA or DMND (diamond binary) format in a specific directory that must be referenced in the mosca command.

```
docker run -v /path/to/folder_of_input/:/input_data -v /path/to/output_folder/:/MOSCA_analysis -v /path/to/folder_of_databases:/MOSCA/Databases/annotation_databases mosca [arguments]
```

The command is to be inputed as presented here, with some alterations:
* "/path/to/folder_of_input/" must be replaced with the **absolute path** to the folder containing the files to be analysed by MOSCA
* "/path/to/output_folder/" must be replaced with the **absolute path** to the folder where to output results from analysis with MOSCA
* "/path/to/folder_of_databases/" must be replaced with the **absolute path** to the folder containing the databases for annotation with DIAMOND
* [arguments] is to be replaced with the parameters exactly as would be inputed when using MOSCA outside of docker. **Filenames (option '-f'), however, must be inputed as /input_data/name_of_file!**

An example command for **Two single-end files (two different experiments)**

```
docker run -v /path/to/my_data_folder/:/input_data -v /path/to/output_folder/:/MOSCA_analysis -v /path/to/my_data_folder/:/input_data -v /path/to/my_database_folder:/MOSCA/Databases/annotation_databases mosca --files /input_data/mg1_R1.fastq,/input_data/mg1_R2.fastq:/input_data/mt1_R1.fastq,/input_data/mt1_R2.fastq [non-mandatory arguments]
```