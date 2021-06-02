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
&& bash MOSCA/workflow/envs/install.bash \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/mosca.py" ]