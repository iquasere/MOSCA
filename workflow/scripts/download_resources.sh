#!/bin/bash

svn export https://github.com/timflutre/trimmomatic/trunk/adapters "$1/adapters"
svn export https://github.com/biocore/sortmerna/trunk/data/rRNA_databases "$1/rRNA_databases"
find "$1/rRNA_databases/*" | grep -v ".fasta" | xargs rm -fr