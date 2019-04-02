FROM continuumio/miniconda:latest

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& cd \home \
&& git clone https://github.com/iquasere/MOSCA.git \
&& wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0.tar.gz \
&& tar -xzf SPAdes-3.9.0.tar.gz \
&& rm SPAdes-3.9.0.tar.gz \
&& export PATH=/home/SPAdes-3.9.0:$PATH \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& conda install fastqc \
&& conda install -c biocore sortmerna \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh \
&& conda install seqtk \
&& conda install -c faircloth-lab trimmomatic \
&& conda install megahit \
&& conda install quast \
&& conda install fraggenescan \
&& conda install diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
&& git clone -b devel https://github.com/claczny/VizBin.git \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda bioconductor-deseq2 \
&& conda install -c bioconda bioconductor-genomeinfodbdata \
&& conda install -c bioconda bioconductor-edger \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install scikit-learn \
&& apt-get purge -y --auto-remove $buildDeps

CMD ['/bin/sh']
