conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba
mamba install bioconductor-deseq2 bioconductor-edger bioconductor-pcamethods bioconductor-rots bioconductor-vsn blast \
bowtie2 checkm-genome diamond fastqc fraggenescan htseq krona maxbin2 megahit openpyxl r-optparse r-pheatmap reportlab \
seqkit seqtk sortmerna=2.1b spades svn trimmomatic upimapi xlrd r-rcolorbrewer pandas scikit-learn lxml biopython \
progressbar33 tqdm xlsxwriter recognizer maxquant quast keggcharter samtools snakemake metaphlan searchgui=3.3 \
peptide-shaker=1.16 python=3 tbb=2020.2 gmcloser -c conda-forge -c bioconda -c anaconda
conda info --base | bash MOSCA/workflow/envs/ci_build.sh
