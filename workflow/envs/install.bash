#!/usr/bin/env bash

conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
mamba env create --file MOSCA/workflow/envs/environment.yml
BASE_DIR="conda info --base"
ENV_DIR="${BASE_DIR}/envs/mosca"
echo "Storing MOSCA's files in the Conda environment at: ${ENV_DIR}"
mkdir -p "${ENV_DIR}/share/MOSCA/scripts" "${ENV_DIR}/bin" "${ENV_DIR}/share/MOSCA/resources"
cp MOSCA/workflow/scripts/* MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${ENV_DIR}/share/MOSCA/scripts"
cp MOSCA/resources/* "${ENV_DIR}/share/MOSCA/resources"
chmod +x "${ENV_DIR}/share/MOSCA/scripts/mosca.py"
ln -s "${ENV_DIR}/share/MOSCA/scripts/mosca.py" "${ENV_DIR}/bin/"
