FROM continuumio/miniconda3:4.7.12

RUN buildDeps='build-essential zlib1g-dev' \
&& mosca_path='/share/MOSCA' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
#&& conda config --add channels r \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git \
&& bash MOSCA/workflow/envs/install.bash
RUN bash MOSCA/workflow/envs/ci_build.sh --conda_path=. --mosca_path=${mosca_path} \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps \
&& bash MOSCA/workflow/envs/ci_build.sh

ENTRYPOINT [ "python", "${mosca_path}/scripts/mosca.py" ]