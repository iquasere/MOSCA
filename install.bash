#-----------FastQC-----------
install https://sourceforge.net/projects/xming/?source=typ_redirect
export DISPLAY=:0
sudo apt-get install -y fastqc
wget https://launchpad.net/ubuntu/+archive/primary/+files/fastqc_0.11.5+dfsg-3_all.deb
sudo dpkg -i fastqc_0.11.5+dfsg-3_all.deb
sudo apt-get install -f

#-----------Anaconda-----------
wget http://repo.continuum.io/archive/Anaconda3-4.0.0-Linux-x86_64.sh
bash Anaconda3-4.0.0-Linux-x86_64.sh
export PATH=~/anaconda3/bin:$PATH

#-----------BMTagger-----------
#conda install -c bioconda bmtagger

#-----------SortMeRNA-----------
conda install -c biocore sortmerna
sudo apt-get install -y seqtk

#-----------Trimmomatic-----------
conda install -c faircloth-lab trimmomatic

sudo apt-get install -y subversion #----> for downloading trimmomatic adapter files

#-----------SPAdes-----------
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0.tar.gz
tar -xzf SPAdes-3.9.0.tar.gz
cd SPAdes-3.9.0
export PATH=/home/jsequeira/SPAdes-3.9.0:$PATH ???
cd ..

#-----------MEGAHIT-----------
sudo apt install aptitude
sudo aptitude install build-essential
sudo apt-get install -y zlib1g-dev
git clone https://github.com/voutcn/megahit.git
cd megahit
make
cd ..

#-----------MetaQUAST-----------
wget https://downloads.sourceforge.net/project/quast/quast-4.5.tar.gz
tar -xzf quast-4.5.tar.gz
cd quast-4.5
./setup.py install_full
cd ..

#-----------FragGeneScan-----------
git clone https://github.com/wltrimbl/FGS.git
cd FGS
make
cd ..

#-----------DIAMOND-----------
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo apt-get install -y python-configobj

#-----------UNIPROT-----------
cd Databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip *
cat *.fasta > uniprot.fasta
cd ..

#-----------VizBin-----------
git clone -b devel https://github.com/claczny/VizBin.git

#-----------Python packages-----------
pip install progressbar33
conda install -c conda-forge biopython
conda install numpy

#-----------R packages-----------
echo -e 'if (!requireNamespace("BiocManager", quietly = TRUE))\n\tinstall.packages("BiocManager")\nBiocManager::install("DESeq2", version = "3.8")' > check_deseq.r
Rscript check_deseq.r

