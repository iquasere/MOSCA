#!/usr/bin/env bash

read -r BASE_DIR
ENV_DIR="${BASE_DIR}/envs/mosca"

echo "Storing MOSCA's files in the Conda environment at: ${ENV_DIR}"
mkdir -p "${ENV_DIR}/share/MOSCA/scripts" "${ENV_DIR}/bin" "${ENV_DIR}/share/MOSCA/resources"
cp MOSCA/workflow/scripts/* MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${ENV_DIR}/share/MOSCA/scripts"
cp -r MOSCA/resources "${ENV_DIR}/share/MOSCA/resources"
chmod +x "${ENV_DIR}/share/MOSCA/scripts/mosca.py"
ln -s "${ENV_DIR}/share/MOSCA/scripts/mosca.py" "${ENV_DIR}/bin/"
