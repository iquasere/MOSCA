FROM continuumio/miniconda3:4.7.12
# should also run with the next version of miniconda image, try that next

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
#&& conda config --add channels r \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git \
&& bash MOSCA/workflow/envs/install.bash \
&& bash MOSCA/workflow/envs/ci_build.sh --conda_path=. --mosca_path=/share/MOSCA \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

ENTRYPOINT [ "python", "bin/mosca.py" ]