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
&& svn export https://github.com/biocore/sortmerna/trunk/data/rRNA_databases /MOSCA/Databases/rRNA_databases \
&& find /MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr \
&& conda install -c bioconda seqtk \
&& conda install -c bioconda trimmomatic \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters /MOSCA/Databases/illumina_adapters \
&& conda install -c bioconda megahit \
&& conda install -c bioconda spades \
# && conda install -c bioconda quast \                                           # TODO - introduce version control so quast can be installed through conda, or wait until it gets python3
&& pip install quast \
&& conda install -c bioconda fraggenescan \
&& conda install -c bioconda diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda checkm-genome \
&& mkdir MOSCA/Databases/checkm_data \
&& cd MOSCA/Databases/checkm_data \
&& curl -L -O https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
&& tar xzf checkm_data_2015_01_16.tar.gz \
&& checkm data setRoot \
&& cd ../../.. \
&& conda install -c anaconda reportlab \
&& conda install -c bioconda bioconductor-deseq2 \
&& conda install -c r r-stringi \
&& conda install -c anaconda openpyxl \
&& conda install -c bioconda bioconductor-edger \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install -c anaconda scikit-learn \
&& conda install -c anaconda lxml \
&& conda install -c bioconda blast \
&& mkdir -p /MOSCA/Databases/annotation_databases \
&& mkdir /input_data \
&& mkdir /output \
&& mkdir /MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz -P /MOSCA/Databases/COG \
&& cd MOSCA/Databases/COG \
&& tar -xzvf cdd.tar.gz --wildcards --no-anchored 'COG*.smp' \
&& cd ../../.. \
&& rm /MOSCA/Databases/COG/cdd.tar.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz -P /MOSCA/Databases/COG \
&& gunzip /MOSCA/Databases/COG/cddid.tbl.gz \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt -P MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog -P MOSCA/Databases/COG \
&& wget https://github.com/aleimba/bac-genomics-scripts/raw/master/cdd2cog/cdd2cog.pl -P MOSCA/scripts \
# sed -i '302s#.*#    my $pssm_id = $1 if $line[1] =~ /^gnl\\|CDD\\|(\\d+)/; \# get PSSM-Id from the subject hit#' MOSCA/cdd2cog.pl      # Sometimes this is needed... will save when use of uninitialized value floods the screen - when they change from CDD:number to gnl|CDD|number
# && conda install -c bioconda searchgui \
# && conda install -c bioconda peptide-shaker \
&& apt-get install -y libpwiz-tools \
&& wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& tar -xzvf SearchGUI-3.3.16-mac_and_linux.tar.gz \
&& wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip \
&& unzip PeptideShaker-1.16.41.zip \
&& conda install -c bioconda -c conda-forge maxquant \
&& conda install -c conda-forge biopython \
&& conda install -c anaconda perl \
&& git clone https://github.com/marbl/Krona.git \
&& cd Krona/KronaTools/ \
&& perl install.pl \
&& apt-get install poppler-utils \
# && conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/MOSCA/scripts/mosca.py" ]