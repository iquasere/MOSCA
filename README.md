![ScreenShot](mosca_logo.png)

Logo by [Sérgio A. Silva](https://www.ceb.uminho.pt/People/Details/64888072-5cde-42da-b7e5-691d380cefb2)

# Meta-Omics Software for Community Analysis (MOSCA)

**MOSCA** (portuguese for fly) is a pipeline designed for performing metagenomics (MG) and metatranscriptomics (MT) integrated data analyses, in a mostly local and fully automated workflow.

## Table of contents

* [Features](#features)
* [Installation](#setting-up-mosca)
  * [with Bioconda](#installation-with-bioconda)
  * [with Docker](#mosca-is-also-available-as-a-docker-image!)
* [Running MOSCA](#base-arguments-for-running-mosca)
* [MOSCA's GUI](#mosguito)

## Features

* **Preprocessing** where low quality regions of data are trimmed and reads less interest are removed. FastQC's reports are used to automatically set the parameters for the other tools. It includes:
    * initial quality check with **FastQC**
    * Illumina artificial sequences removal with **Trimmomatic**: based on **FastQC** reports, MOSCA will find the adapters file most approprita to the data
    * rRNA removal with **SortMeRNA**: uses Pfam and SILVA databases as reference
    * quality trimming with **Trimmomatic**: 
        * another **FastQC** report will be generated after rRNA removal, and will be used to set the parameters for **Trimmomatic**'s hard trimmers (CROP and HEADCROP). This will ensure that the data will be reported as excellent by FastQC
        * reads with less than 20 average quality or 100 nuleotides of length will also be removed
    * final quality check with **FastQC** 
* **Assembly** where MG trimmed reads will be assembled to partially reconstruct the original genomes in the samples. It includes:
    * assembly with two possible assemblers - **MetaSPAdes** and **Megahit** - which will be used in a multi-kmer approach
    * control over the quality of the contigs, with **MetaQUAST** reporting on several classical metrics (such as N50 and L50) and alignment of reads for estimating percentage of reads used in assembly, with **Bowtie2**
* **Annotation** where proteins present in the contigs will be identified. It includes:
    * gene calling with **FragGeneScan**
    * annotation of identified ORFs with **DIAMOND**, using the **UniProt database** as reference - MOSCA only reports on the first annotation
    * retrieval of diverse biological information with [**UPIMAPI**](https://anaconda.org/bioconda/UPIMAPI)
    * functional annotation with [**reCOGnizer**](https://anaconda.org/bioconda/reCOGnizer), using the **COG database** as reference
        * MOSCA automatically **generates new databases by the number of threads specified**, thus allowing for multithread annotation with **RPSBLAST**
    * the quantification of each protein in MG data, by alignment of MG reads to the contigs using **Bowtie2** and quantification of reads to protein using **HTSeq-count**
* **Binning** where the contigs are clustered into taxonomic units, to validate (or not) the annotation, and possibily help reconstructing genomes from the samples
    * **MaxBin2** bins the contigs by tetranucleotide composition, relative abundance, and marker genes analysis
    * the final bins are reported for their completeness - how many of the marker genes are present in each bin
* **Metatranscriptomics analysis** where the expression of each identified gene is quantified. It includes:
    * alignment of MT reads to the MG contigs with **Bowtie2**, and quantification of reads to protein using **HTSeq-count**
    * differential gene expression and multisample comparison using **DESeq2**
* **Normalization** of protein quantification for the final reports using **edgeR**
* **Pathway representation** with [**KEGGCharter**](https://anaconda.org/bioconda/KEGGCharter), representing both the metabolic networks of most abundant taxa and expression levels of metabolic functions

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

### Installation with Bioconda

To use MOSCA through Bioconda, an environment must be created containing all the dependencies. 
Create the environment and install MOSCA with 
```
conda create -n mosca -c conda-forge -c bioconda -c anaconda mosca
```
Activate the environment with 
```
conda activate mosca
```
If that command is over without error, you have successfully installed MOSCA!

### MOSCA is also available as a Docker image!

To use MOSCA's Docker version, Docker must first be installed.

```
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```

After having docker installed in the system, pull MOSCA's image.

```
docker pull iquasere/mosca:latest
```

## Base arguments for running MOSCA

Since 1.2.0, MOSCA only accepts input from a config file, in either JSON or YAML format.
This repo has an available [config file](https://github.com/iquasere/MOSCA/blob/development/config/config.json), 
which can be used for MOSCA as follows:
```
python mosca.py --configfile config.json
```
The config file allows to customize MOSCA's workflow, but for the convenience of users, many typical decisions in MG and
MT workflow are already automized. The customization, therefore, is only related to steps that are not yet well established
in the field of MG (e.g. assembling data into contigs is still a controversial step that may lose information on data).

Following are the options available in the config file, and the accepted values:

|            Parameter           |                            Options                           | Required | Description                                                                                                                                                                                                                   |
|:------------------------------:|:------------------------------------------------------------:|:--------:|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|             output             |                            String                            |    Yes   | Name of folder where MOSCA's results will be stored (if it doesn't exist, it will be created)                                                                                                                                 |
|             threads            |                              Int                             |    Yes   | Number of maximum threads for MOSCA to use                                                                                                                                                                                    |
|           experiments          |                            String                            |    Yes   | Name of TSV file with information on samples/files/conditions                                                                                                                                                                 |
| trimmomatic_adapters_directory |                            String                            |    Yes   | Name of folder containing adapters for Trimmomatic's ADAPTER REMOVAL preprocessing tool                                                                                                                                       |
|    rrna_databases_directory    |                            String                            |    Yes   | Name of folder containing rRNA databases to use as reference for rRNA removal with SortMeRNA                                                                                                                                  |
|            assembler           |                      metaspades, megahit                     |    Yes   | Name of assembler to use for iterative co-assembly of MG data                                                                                                                                                                 |
|            markerset           |                            40, 107                           |    Yes   | Name of markerset to use for completeness/contamination estimation with CheckM over the contigs obtained with MaxBin2                                                                                                         |
|           error_model          | sanger_5, sanger_10, 454_10, 454_30, illumina_5, illumina_10 |    No    | Name of file to use as the error model for gene calling with FragGeneScan. sanger, 454 or illumina if either Sanger, pyro- or Illumina sequencing reads are the input to gene calling. Leave empty if assembly was performed. |
|        diamond_database        |                            String                            |    Yes   | Name of FASTA or DMND (DIAMOND formatted database) file to use as input for annotation with DIAMOND                                                                                                                           |
|        download_uniprot        |                          TRUE, FALSE                         |    Yes   | If UniProtKB (SwissProt + TrEMBL) is to be download. If TRUE, will download it to the folder indicated in diamond_database                                                                                                    |
|     diamond_max_target_seqs    |                              Int                             |    Yes   | Number of matches to report for each protein from annotation with DIAMOND                                                                                                                                                     |
| recognizer_databases_directory |                            String                            |    Yes   | Name of folder containing the resources for reCOGnizer annotation. If those are not present in the folder, they will be downloaded                                                                                            |
|      normalization_method      |                           TMM, RLE                           |    Yes   | Method to use for normalization                                                                                                                                                                                               |
|        keggcharter_maps        |            Comma-separated list of KEGG maps' IDs            |    No    | If empty, KEGGCharter will use the default prokaryotic maps. These metabolic maps will have MG information represented in them, and gene expression if MT data is available                                                   |
|     keggcharter_taxa_level     |  SPECIES, GENUS, FAMILY, ORDER, CLASS, PHYLUM, SUPERKINGDOM  |    Yes   | The taxonomic level to represent with KEGGCharter. If above SPECIES, KEGGCharter will represent group information and represent is as such for each taxonomic group                                                           |
|   keggcharter_number_of_taxa   |                     Int, ideally under 11                    |    Yes   | How many of the most abundant taxa should be represented with KEGGCharter                                                                                                                                                     |
|    reporter_lists_directory    |                            String                            |    Yes   | Name of folder containing lists for reporter module of MOSCA                                                                                                                                                                  |

## MOSGUITO

MOSca's GUI TO generate config files ([MOSGUITO](https://iquasere.github.io/MOSGUITO/)) can be used to produce a configuration file for MOSCA. In it you can also find the [link](https://docs.google.com/spreadsheets/d/12BppOf32QPRl6Ey39ACWoOXPVWtTrGIUN2u-_BpwSvo/edit#gid=0) to an experiment's file, required to input information for each data file to be analysed with MOSCA. 