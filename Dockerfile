FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git \
&& conda install svn reportlab openpyxl xlrd>=0.9.0 r-rcolorbrewer pandas scikit-learn lxml biopython perl \
&& conda install -c bioconda fastqc sortmerna=2.1 seqtk trimmomatic megahit spades fraggenescan diamond upimapi htseq bowtie2 maxbin2 checkm-genome bioconductor-deseq2=1.22.1 bioconductor-edger=3.24.3 r-pheatmap r-optparse blast krona seqkit samtools \
&& conda install -c conda-forge progressbar33 tqdm>=4.33.0 xlsxwriter unzip \
&& conda install -c bioconda -c conda-forge recognizer maxquant keggcharter quast \
&& conda clean --all \
&& perl /opt/conda/opt/krona/install.pl \
&& wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip \
&& unzip PeptideShaker-1.16.41.zip \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/Databases/illumina_adapters \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/MOSCA/scripts/mosca.py" ]