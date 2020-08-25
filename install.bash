#!/usr/bin/env bash

buildDeps='build-essential zlib1g-dev'
apt-get update
apt-get install -y $buildDeps --no-install-recommends

# Tools installation
source ~/anaconda3/etc/profile.d/conda.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y svn reportlab openpyxl xlrd>=0.9.0 r-rcolorbrewer pandas scikit-learn lxml biopython perl
conda install -y -c bioconda fastqc sortmerna=2.1 seqtk trimmomatic megahit spades fraggenescan diamond upimapi htseq bowtie2 checkm-genome bioconductor-edger=3.24.3 r-pheatmap r-optparse blast krona seqkit samtools
conda install -y -c bioconda maxbin2 
conda install -y -c bioconductor-deseq2=1.22.1
conda install -y -c conda-forge progressbar33 tqdm>=4.33.0 xlsxwriter unzip
conda install -y -c bioconda -c conda-forge recognizer maxquant quast keggcharter
apt-get install -y libpwiz-tools poppler-utils
perl ~/anaconda3/opt/krona/install.pl
wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz
tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz
wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip
unzip PeptideShaker-1.16.41.zip

# Databases download
svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/Databases/illumina_adapters
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > MOSCA/Databases/annotation_databases/uniprot.fasta
rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
conda clean --all