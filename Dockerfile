FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git -b development \
&& conda install -c bioconda fastqc \
&& conda install -c biocore sortmerna=2.1 \
&& conda install -c anaconda svn \
&& conda install -c bioconda seqtk \
&& conda install -c bioconda trimmomatic \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters /MOSCA/Databases/illumina_adapters \
&& conda install -c bioconda megahit \
&& conda install -c bioconda spades \
&& pip install quast \
&& conda install -c bioconda fraggenescan \
&& conda install -c bioconda diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda checkm-genome \
&& conda install -c anaconda reportlab \
&& conda install -c bioconda bioconductor-deseq2=1.22.1 \
&& conda install -c anaconda openpyxl \
&& conda install -c bioconda bioconductor-edger=3.24.3 \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install -c anaconda scikit-learn \
&& conda install -c anaconda lxml \
&& conda install -c bioconda blast \
&& mkdir -p /MOSCA/Databases/annotation_databases \
&& apt-get update \
&& apt-get install -y libpwiz-tools \
&& wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& apt-get install -y zip unzip \
&& wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip \
&& unzip PeptideShaker-1.16.41.zip \
&& conda install -c bioconda -c conda-forge maxquant \
&& conda install -c conda-forge biopython \
&& conda install -c anaconda perl \
&& git clone https://github.com/marbl/Krona.git \
&& apt-get install -y poppler-utils \
&& git clone https://github.com/iquasere/UPIMAPI.git \
&& git clone https://github.com/iquasere/reCOGnizer.git \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P reCOGnizer \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/MOSCA/scripts/mosca.py" ]