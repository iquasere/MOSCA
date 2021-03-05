FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& conda config --add channels defaults \
&& conda config --add channels bioconda \
&& conda config --add channels conda-forge \
&& git clone https://github.com/iquasere/MOSCA.git \
# obtained with conda env export --no-builds > environment.yml
&& conda env create -f MOSCA/workflow/envs/environment.yml \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps \
&& dir="/share/MOSCA" \
&& mkdir -p "${dir}/scripts" "/bin" \
&& cp MOSCA/workflow/scripts/* "${dir}/scripts" \
&& cp MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${dir}/scripts" \
&& cp -r resources "${dir}/resources" \
&& chmod +x "${dir}/scripts/mosca.py" \
&& ln -s "${dir}/scripts/mosca.py" "${PREFIX}/bin/"

ENTRYPOINT [ "python", "${dir}/scripts/mosca.py" ]