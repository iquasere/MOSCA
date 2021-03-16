#!/bin/bash

svn export https://github.com/timflutre/trimmomatic/trunk/adapters "$1/adapters" --force
svn export https://github.com/biocore/sortmerna/trunk/data/rRNA_databases "$1/rRNA_databases" --force