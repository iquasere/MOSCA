FROM continuumio/miniconda3:latest

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git \
&& conda install fastqc \
&& conda install -c biocore sortmerna \
&& conda install -c anaconda svn \
&& svn export https://github.com/biocore/sortmerna/trunk/rRNA_databases /MOSCA/Databases/rRNA_databases \
&& find /MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/merge-paired-reads.sh -P /MOSCA/scripts \
&& wget https://github.com/biocore/sortmerna/raw/master/scripts/unmerge-paired-reads.sh -P /MOSCA/scripts \
&& conda install seqtk \
&& conda install -c bioconda trimmomatic \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters /MOSCA/Databases/illumina_adapters \
&& conda install megahit \
&& conda install -c bioconda spades \
&& conda install quast \
&& conda install fraggenescan \
&& conda install diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
#&& git clone -b devel https://github.com/claczny/VizBin.git \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda bioconductor-deseq2 \
&& conda install -c bioconda bioconductor-genomeinfodbdata \
&& conda install bioconductor-edger \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install scikit-learn \
&& conda install -y -c bioconda blast \
&& mkdir /MOSCA/Databases/annotation_databases \
&& mkdir -p /MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz -P /MOSCA/Databases/COG \
&& tar -xzvf /MOSCA/Databases/COG/cdd.tar.gz --wildcards --no-anchored 'COG*.smp' -C /MOSCA/Databases/COG \
&& rm /MOSCA/Databases/COG/cdd.tar.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P /MOSCA/Databases/COG \
&& gunzip /MOSCA/Databases/COG/cddid.tbl.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P MOSCA/Databases/COG \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P MOSCA/scripts \
#sed -i '302s#.*#    my $pssm_id = $1 if $line[1] =~ /^gnl\\|CDD\\|(\\d+)/; \# get PSSM-Id from the subject hit#' MOSCA/cdd2cog.pl      # Sometimes this is needed... will save when use of uninitialized value floods the screen - when they change from CDD:number to gnl|CDD|number
&& conda install -c bioconda searchgui \
&& conda install -c bioconda peptide-shaker \
&& conda install -c bioconda maxquant \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

ENV files None
ENV data paired
ENV assembler metaspades
ENV annotation_database Databases/annotation_databases/uniprot.fasta
ENV output MOSCA_analysis
ENV output_level max
ENV conditions None
ENV threads 1
ENV memory None
ENV marker_gene_set 40

ENTRYPOINT [ "python", "/MOSCA/mosca.py" ]