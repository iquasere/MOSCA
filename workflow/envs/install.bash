conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba
mamba env create --file MOSCA/workflow/envs/environment.yml
conda info --base | bash MOSCA/workflow/envs/ci_build.sh -
