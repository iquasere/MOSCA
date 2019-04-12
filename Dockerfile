FROM continuumio/miniconda:latest

EXPOSE 5000

WORKDIR /MOSCA



RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& cd /home \
&& git clone https://github.com/iquasere/MOSCA.git \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& conda install fastqc \
&& conda install -c biocore sortmerna \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh -P MOSCA \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh -P MOSCA \
&& conda install -c anaconda svn \
&& mkdir -p MOSCA/Databases/rRNA_databases \
&& svn checkout https://github.com/biocore/sortmerna/trunk/rRNA_databases MOSCA/Databases/rRNA_databases \
&& find MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr \
&& conda install seqtk \
&& conda install -c bioconda trimmomatic \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters MOSCA/Databases/illumina_adapters \
&& conda install megahit \
&& conda install -c bioconda spades \
&& conda install quast \
&& conda install fraggenescan \
&& conda install diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
&& git clone -b devel https://github.com/claczny/VizBin.git \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda bioconductor-deseq2 \
&& conda install -c bioconda bioconductor-edger \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install scikit-learn \
&& conda install -c anaconda lxml \
&& conda install -c bioconda searchgui \
&& conda install -c bioconda peptide-shaker \
&& conda install -c bioconda maxquant \
&& mkdir -p MOSCA/Databases/COG \
# COGs involve over 300Mb of data, should be included?
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P MOSCA/Databases/COG \
&& gunzip MOSCA/Databases/COG/cddid.tbl.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz -P MOSCA/Databases/COG \
&& tar -xvzf MOSCA/Databases/COG/Cog_LE.tar.gz -C MOSCA/Databases/COG \
&& rm MOSCA/Databases/COG/Cog_LE.tar.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P MOSCA/Databases/COG \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P MOSCA \
&& sed -i '302s#.*#    my $pssm_id = $1 if $line[1] =~ /^gnl\|CDD\|(\d+)/; \# get PSSM-Id from the subject hit#' MOSCA/cdd2cog.pl \
&& apt-get purge -y --auto-remove $buildDeps

CMD ['/bin/sh']