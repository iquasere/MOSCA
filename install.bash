#!/usr/bin/env bash

# Proteomics tools installation
apt-get install -y libpwiz-tools poppler-utils
perl ~/anaconda3/opt/krona/install.pl
wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz
tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz
wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip
unzip PeptideShaker-1.16.41.zip

# Databases download
svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/resources/illumina_adapters
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > MOSCA/resources/uniprot.fasta
rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz