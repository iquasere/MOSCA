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
&& conda install seqtk \
&& conda install -c bioconda trimmomatic \
&& svn export https://github.com/timflutre/trimmomatic/trunk/adapters /MOSCA/Databases/illumina_adapters \
&& conda install megahit \
&& conda install -c bioconda spades \
# && conda install -c bioconda quast \                                            # TODO - introduce version control so quast can be installed through conda, or wait until it gets python3
&& pip install quast \
&& conda install -c bioconda fraggenescan \
&& conda install -c bioconda diamond \
&& conda install -c conda-forge progressbar33 \
&& conda install -c bioconda htseq \
&& conda install -c bioconda bowtie2 \
&& conda install -c bioconda maxbin2 \
&& conda install -c bioconda checkm \
&& curl -L -O https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz \
&& tar xzf checkm_data_2015_01_16.tar.gz \
&& conda create -n py27 python=2.7 \
&& echo "source activate py27" > ~/.bashrc
#ENV PATH /opt/conda/envs/env/bin:$PATH 											# CheckM installation is still a nono
#CMD [ " conda activate py27 && checkm data setRoot ." ]
RUN conda install -c anaconda biopython \
&& conda install -c anaconda reportlab \
&& conda install -c bioconda bioconductor-deseq2 \
&& conda install -c r r-stringi \                                               # reference to https://github.com/jupyter/docker-stacks/issues/927 (loading DESeq2 fails otherwise)
&& conda install -c anaconda openpyxl \                                         # normalization fails otherwise with "No module named 'openpyxl'"
&& conda install bioconductor-edger \
&& conda install -c bioconda r-pheatmap \
&& conda install -c r r-rcolorbrewer \
&& conda install -c bioconda r-optparse \
&& conda install -c anaconda pandas \
&& conda install -c conda-forge tqdm \
&& conda install scikit-learn \
&& conda install -c bioconda blast \
&& mkdir -p /MOSCA/Databases/annotation_databases \
&& mkdir /input_data \
&& mkdir /output \
&& mkdir /MOSCA/Databases/COG \
&& wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz -P /MOSCA/Databases/COG \
# && tar -xzvf /MOSCA/Databases/COG/cdd.tar.gz --wildcards --no-anchored 'COG*.smp' -C /MOSCA/Databases/COG \
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
&& conda install -c anaconda lxml \
&& conda install -c bioconda searchgui \
&& conda install -c bioconda peptide-shaker \
&& conda install -c bioconda -c conda-forge maxquant \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "/MOSCA/scripts/mosca.py" ]