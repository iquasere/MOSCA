buildDeps='build-essential zlib1g-dev'
apt-get update
apt-get install -y $buildDeps --no-install-recommends
rm -rf /var/lib/apt/lists/*
cd ~
git clone https://github.com/iquasere/MOSCA.git
wget http://spades.bioinf.spbau.ru/release3.11.1/SPAdes-3.11.1.tar.gz
tar -xzf SPAdes-3.11.1.tar.gz
rm SPAdes-3.11.1.tar.gz
export PATH=/home/SPAdes-3.9.0:$PATH
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y fastqc
conda install -y -c biocore sortmerna
wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh
wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh
conda install -y seqtk
conda install -y -c faircloth-lab trimmomatic
conda install -y megahit
conda install -y quast
conda install -y fraggenescan
conda install -y diamond
conda install -y -c conda-forge progressbar33
conda install -y -c bioconda htseq
conda install -y -c bioconda bowtie2
git clone -b devel https://github.com/claczny/VizBin.git
conda install -y -c bioconda maxbin2
conda install -y -c bioconda bioconductor-deseq2
R -e 'BiocManager::install("GenomeInfoDbData", version = "3.8")'
conda install -y -c bioconda bioconductor-genomeinfodbdata
conda install -y -c bioconda bioconductor-edger
conda install -y -c bioconda r-pheatmap
conda install -y -c r r-rcolorbrewer
conda install -y -c bioconda r-optparse
conda install -y -c anaconda pandas
conda install -y -c conda-forge tqdm
conda install -y scikit-learn
mkdir MOSCA/databases
cd MOSCA/databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > uniprot.fasta.gz
rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
gunzip uniprot.fasta.gz
