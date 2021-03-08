#!/usr/bin/env bash

while [ $# -gt 0 ]; do
  case "$1" in
    --conda_path=*)
      conda_path="${1#*=}"
      ;;
    --mosca_path=*)
      mosca_path="${1#*=}"
      ;;
    *)
      printf "***************************\n"
      printf "* Error: Invalid argument.*\n"
      printf "***************************\n"
      exit 1
  esac
  shift
done

mkdir -p "${mosca_path}/scripts" "${conda_path}/bin"
cp MOSCA/workflow/scripts/* "${mosca_path}/scripts"
cp MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py "${mosca_path}/scripts"
cp -r MOSCA/resources "${mosca_path}/resources"
chmod +x "${mosca_path}/scripts/mosca.py"
rm "${conda_path}/bin/mosca.py"
ln -s "${mosca_path}/scripts/mosca.py" "${conda_path}/bin/"
