buildDeps='build-essential zlib1g-dev'
apt-get update
apt-get install -y $buildDeps zlib1g-dev --no-install-recommends
-rf /var/lib/apt/lists/*
cd ~
git clone https://github.com/iquasere/MOSCA.git
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0.tar.gz
tar -xzf SPAdes-3.9.0.tar.gz
rm SPAdes-3.9.0.tar.gz
export PATH=/home/SPAdes-3.9.0:$PATH
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install fastqc
conda install -c biocore sortmerna
wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh -P MOSCA
wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh -P MOSCA
conda install seqtk
conda install -c faircloth-lab trimmomatic
conda install megahit
conda install quast
conda install fraggenescan
conda install diamond
conda install -c conda-forge progressbar33
conda install -c bioconda htseq
conda install -c bioconda bowtie2
conda install -c r r
conda install -c bioconda bioconductor-deseq2
conda install -c bioconda bioconductor-edger
conda install -c bioconda r-pheatmap
Rscript /MOSCA/install_r_packages.R
#docker pull 990210oliver/mycc.docker:v1                                        not working... yet!
git clone -b devel https://github.com/claczny/VizBin.git
mkdir MOSCA/databases
cd MOSCA/databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > uniprot.fasta.gz
rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
gunzip uniprot.fasta.gz
