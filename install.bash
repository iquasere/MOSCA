#!/bin/bash

# Set default values for conda_dir and mosca_dir
conda_dir=$(conda info --base)
mosca_env="${conda_dir}/envs/mosca"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --conda_dir)
      conda_dir="$2"
      shift # past argument
      shift # past value
      ;;
    --mosca_env)
      mosca_env="$2"
      shift # past argument
      shift # past value
      ;;
  esac
done

# Set default value for mosca_dir if conda_dir is set but mosca_dir is not
if [ -z "mosca_env" ] && [ ! -z "$conda_dir" ]; then
  mosca_env="${conda_dir}/envs/mosca"
fi

conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge

echo "Storing MOSCA's files in the Conda environment at: ${mosca_env}"
mkdir -p "${mosca_env}/share/MOSCA/scripts" "${mosca_env}/bin" "${mosca_env}/share/MOSCA/resources"
cp MOSCA/workflow/scripts/* MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${mosca_env}/share/MOSCA/scripts"
cp MOSCA/resources/* "${mosca_env}/share/MOSCA/resources"
chmod +x "${mosca_env}/share/MOSCA/scripts/mosca.py"
ln -s "${mosca_env}/share/MOSCA/scripts/mosca.py" "${mosca_env}/bin/"