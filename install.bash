#!/usr/bin/env bash

dir="${PREFIX}/share/MOSCA/resources"

# Databases download
svn export https://github.com/timflutre/trimmomatic/trunk/adapters "${dir}/illumina_adapters"
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz > "${dir}/uniprot.fasta"
rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz
svn export https://github.com/biocore/sortmerna/trunk/rRNA_databases "${dir}/rRNA_databases"
find "${dir}/rRNA_databases/*" | grep -v ".fasta" | xargs rm -fr

# Proteomics tools installation
apt-get install -y libpwiz-tools poppler-utils
perl ~/anaconda3/opt/krona/install.pl
wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz
tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz
wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip
unzip PeptideShaker-1.16.41.zip