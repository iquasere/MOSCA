FROM continuumio/miniconda3:4.9.2
# should also run with the next version of miniconda image, try that next

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& git clone https://github.com/iquasere/MOSCA.git \
&& conda config --add channels anaconda \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
#&& conda env create -f MOSCA/workflow/envs/environment.yml \
&& conda install -c conda-forge mamba \
&& mamba install bioconductor-deseq2 bioconductor-edger bioconductor-pcamethods bioconductor-rots bioconductor-vsn blast bowtie2 checkm-genome diamond fastqc fraggenescan htseq krona maxbin2 megahit openpyxl r-optparse r-pheatmap reportlab seqkit seqtk sortmerna spades svn trimmomatic upimapi xlrd r-rcolorbrewer pandas scikit-learn lxml biopython progressbar33 tqdm xlsxwriter recognizer maxquant quast keggcharter samtools snakemake metaphlan searchgui peptide-shaker python=3 -c conda-forge -c bioconda -c anaconda \
&& bash MOSCA/workflow/envs/ci_build.sh --conda_path=. --mosca_path=/share/MOSCA \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/mosca.py" ]