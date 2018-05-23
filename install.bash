#-----------FastQC-----------
install https://sourceforge.net/projects/xming/?source=typ_redirect
export DISPLAY=:0
sudo apt-get install fastqc
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
wget ./merge-paired-reads.sh ./unmerge-paired-reads.sh !!!!
sed -i 's/\r$//' ./merge-paired-reads.sh

#-----------Trimmomatic-----------
conda install -c faircloth-lab trimmomatic
PATH=$PATH:/home/sequeira/anaconda3/jar

sudo apt-get install subversion ----> for downloading trimmomatic adapter files

#-----------SPAdes-----------
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0.tar.gz
tar -xzf SPAdes-3.9.0.tar.gz
cd SPAdes-3.9.0
export PATH=/home/jsequeira/SPAdes-3.9.0:$PATH ???

#-----------MEGAHIT-----------
sudo apt install aptitude
sudo aptitude install build-essential
sudo apt-get install zlib1g-dev
git clone https://github.com/voutcn/megahit.git
cd megahit
make

#-----------MetaQUAST-----------
wget https://downloads.sourceforge.net/project/quast/quast-4.5.tar.gz
tar -xzf quast-4.5.tar.gz
cd quast-4.5
./setup.py install_full

#-----------FragGeneScan-----------
git clone https://github.com/wltrimbl/FGS.git
cd FGS
make

#-----------DIAMOND-----------
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.9/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo apt-get install python-configobj

#-----------UNIPROT-----------
cd Databases
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
cat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > uniprot.fasta.gz
cd ..

#-----------Python packages-----------
pip install progressbar33

import numpy as np
import subprocess