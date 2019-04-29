#!/usr/bin/env bash

buildDeps='build-essential zlib1g-dev'
apt-get update
apt-get install -y $buildDeps --no-install-recommends
source ~/anaconda3/etc/profile.d/conda.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y fastqc
conda install -y -c biocore sortmerna
conda install -y -c anaconda svn
svn export https://github.com/biocore/sortmerna/trunk/rRNA_databases MOSCA/Databases/rRNA_databases    # TODO - might have problem, svn not found
find MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr
wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh -P MOSCA
wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh -P MOSCA
conda install -y seqtk
conda install -y -c faircloth-lab trimmomatic
svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/Databases/illumina_adapters    # TODO - might have problem, svn not found
conda install -y megahit
conda install -y -c bioconda spades
conda install -y quast
conda install -y fraggenescan
conda install -y diamond
conda install -y -c conda-forge progressbar33
conda install -y -c bioconda htseq
conda install -y -c bioconda bowtie2
git clone -b devel https://github.com/claczny/VizBin.git
conda install -y -c bioconda maxbin2
conda install -y -c bioconda bioconductor-deseq2=1.22.1
conda install -y -c bioconda bioconductor-genomeinfodbdata=1.16.0
conda install -y -c bioconda bioconductor-edger
conda install -y -c bioconda r-pheatmap
conda install -y -c r r-rcolorbrewer
conda install -y -c bioconda r-optparse
conda install -y -c anaconda pandas
conda install -y -c conda-forge tqdm
conda install -y scikit-learn
conda install -c anaconda lxml
mkdir MOSCA/Databases/annotation_databases
cd MOSCA/Databases/annotation_databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip *.fasta.gz
cat uniprot_trembl.fasta uniprot_sprot.fasta > uniprot.fasta
rm uniprot_trembl.fasta uniprot_sprot.fasta
cd ../../..
mkdir -p MOSCA/Databases/COG
wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P MOSCA/Databases/COG
gunzip MOSCA/Databases/COG/cddid.tbl.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz -P MOSCA/Databases/COG
tar -xvzf MOSCA/Databases/COG/Cog_LE.tar.gz -C MOSCA/Databases/COG
rm MOSCA/Databases/COG/Cog_LE.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P MOSCA/Databases/COG
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P MOSCA/Databases/COG
wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P MOSCA
sed -i '302s#.*#    my $pssm_id = $1 if $line[1] =~ /^gnl\|CDD\|(\d+)/; \# get PSSM-Id from the subject hit#' MOSCA/cdd2cog.pl
conda install -y -c bioconda searchgui
conda install -y -c bioconda peptide-shaker
conda install -y -c bioconda maxquant

