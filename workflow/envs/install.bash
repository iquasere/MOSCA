conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y -c conda-forge mamba
mamba env create --file MOSCA/workflow/envs/environment.yml
which python | awk '{print substr($0, 1, length($0)-6)}' | xargs bash MOSCA/workflow/envs/ci_build.sh --conda_path=. --mosca_path=/share/MOSCA