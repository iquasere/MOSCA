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
  shift # past argument
done

# Set default value for mosca_dir if it is not set
if [ -z "$mosca_env" ]; then
  mosca_env="${conda_dir}/envs/mosca"
fi

conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge

echo "Storing MOSCA's files in the Conda environment at: ${mosca_env}"
# create folders for storing MOSCA's YAMLs and scripts
mkdir -p "${mosca_env}/share/MOSCA" "${mosca_env}/bin"
# copy YAMLs and scripts to the MOSCA Conda environment
cp -r MOSCA/workflow/scripts/* MOSCA/workflow/Snakefile MOSCA/workflow/mosca.py MOSCA/resources MOSCA/workflow/envs  \
  "${mosca_env}/share/MOSCA"
# make MOSCA's main script executable
chmod +x "${mosca_env}/share/MOSCA/scripts/mosca.py"
# create a symbolic link to MOSCA's main script in the bin folder
ln -s "${mosca_env}/share/MOSCA/scripts/mosca.py" "${mosca_env}/bin/"
