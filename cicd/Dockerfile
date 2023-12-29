FROM continuumio/miniconda3

RUN buildDeps='build-essential zlib1g-dev' \
&& apt-get update \
&& apt-get install -y $buildDeps --no-install-recommends \
&& rm -rf /var/lib/apt/lists/* \
&& git clone https://github.com/iquasere/MOSCA.git \
&& bash MOSCA/install.bash \
&& conda clean --all \
&& apt-get purge -y --auto-remove $buildDeps

CMD [ "python", "bin/mosca.py" ]