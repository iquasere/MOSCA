#!/usr/bin/env bash

buildDeps='build-essential zlib1g-dev'
apt-get update
apt-get install -y $buildDeps --no-install-recommends
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c bioconda fastqc
conda install -y -c bioconda sortmerna=2.1
conda install -y -c anaconda svn
svn export https://github.com/biocore/sortmerna/trunk/data/rRNA_databases MOSCA/Databases/rRNA_databases
find MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr
conda install -y -c bioconda seqtk
conda install -y -c bioconda trimmomatic
svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/Databases/illumina_adapters
conda install -y -c bioconda megahit
conda install -y -c bioconda spades
# conda install -y -c bioconda quast                                            # TODO - introduce version control so quast can be installed through conda
pip install quast
conda install -y -c bioconda fraggenescan
conda install -y -c bioconda diamond
conda install -y -c conda-forge progressbar33
conda install -y -c bioconda htseq
conda install -y -c bioconda bowtie2
conda install -y -c bioconda maxbin2
conda install -y -c bioconda checkm-genome
mkdir MOSCA/Databases/checkm_data
cd MOSCA/Databases/checkm_data
curl -L -O https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xzf checkm_data_2015_01_16.tar.gz
checkm data setRoot .
cd ../../..
conda install -y -c anaconda reportlab
conda install -y -c bioconda bioconductor-deseq2=1.22.1
conda install -y -c anaconda openpyxl                                           # normalization fails otherwise with "No module named 'openpyxl'"
#conda install -y -c bioconda bioconductor-genomeinfodb                         # genomeinfodb might be required after deseq
#conda install -y -c bioconda bioconductor-genomeinfodbdata                     # genomeinfodbdata might be required after deseq
conda install -y -c bioconda bioconductor-edger=3.24.3
conda install -y -c bioconda r-pheatmap
conda install -y -c anaconda r-rcolorbrewer
conda install -y -c bioconda r-optparse
conda install -y -c anaconda pandas
conda install -y -c anaconda xlrd                                               # pandas Excel support
conda install -y -c conda-forge tqdm
conda install -y -c anaconda scikit-learn
conda install -y -c anaconda lxml
conda install -y -c bioconda blast
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
cat uniprot_trembl.fasta uniprot_sprot.fasta > MOSCA/Databases/annotation_databases/uniprot.fasta
rm uniprot_trembl.fasta uniprot_sprot.fasta uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
mkdir -p MOSCA/Databases/COG
wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz -P MOSCA/Databases/COG
# && tar -xzvf /MOSCA/Databases/COG/cdd.tar.gz --wildcards --no-anchored 'COG*.smp' -C /MOSCA/Databases/COG
cd MOSCA/Databases/COG
tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp'
cd ../../..
wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P MOSCA/Databases/COG
gunzip MOSCA/Databases/COG/cddid.tbl.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P MOSCA/Databases/COG
wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P MOSCA/Databases/COG
wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P MOSCA/scripts
#sed -i '302s#.*#    my $pssm_id = $1 if $line[1] =~ /^gnl\\|CDD\\|(\\d+)/; \# get PSSM-Id from the subject hit#' MOSCA/cdd2cog.pl      # Sometimes this is needed... will save when use of uninitialized value floods the screen - when they change from CDD:number to gnl|CDD|number
# Proteomics
# conda install -y -c bioconda searchgui
# conda install -y -c bioconda peptide-shaker
apt-get install -y libpwiz-tools
wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz
tar -xzvf SearchGUI-3.3.16-mac_and_linux.tar.gz
wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip
unzip PeptideShaker-1.16.41.zip
conda install -y -c bioconda -c conda-forge maxquant
conda install -y -c anaconda biopython
# Krona plotting
conda install -y -c anaconda perl
git clone https://github.com/marbl/Krona.git
cd Krona/KronaTools/
perl install.pl
apt-get install -y poppler-utils
# conda clean --all