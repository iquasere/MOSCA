#!/usr/bin/env bash

mosca_path="/share/MOSCA"
conda_path="/opt/conda"

mkdir -p "${mosca_path}/scripts" "${conda_path}/bin"
cp MOSCA/workflow/scripts/* "${mosca_path}/scripts"
cp MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${mosca_path}/scripts"
cp MOSCA/resources/* "${mosca_path}/resources"
chmod +x "${mosca_path}/scripts/mosca.py"
ln -s "${mosca_path}/scripts/mosca.py" "${conda_path}/bin/"