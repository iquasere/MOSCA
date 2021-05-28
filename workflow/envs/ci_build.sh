#!/usr/bin/env bash

read -r ENV_DIR

echo "Storing MOSCA's files in the Conda environment at: ${ENV_DIR}"
mkdir -p "${ENV_DIR}/scripts" "${ENV_DIR}/bin"
cp MOSCA/workflow/scripts/* MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${ENV_DIR}/scripts"
cp -r MOSCA/resources "${ENV_DIR}/resources"
chmod +x "${ENV_DIR}/scripts/mosca.py"
ln -s "${ENV_DIR}/scripts/mosca.py" "${ENV_DIR}/bin/"
