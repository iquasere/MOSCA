dir="${PREFIX}/share/MOSCA"
mkdir -p "${dir}/scripts" "${PREFIX}/bin"
cp workflow/scripts/* "${dir}/scripts"
cp workflow/Snakefile workflow/mosca.py "${dir}/scripts"
cp -r resources "${dir}/resources"
chmod +x "${dir}/scripts/mosca.py"
ln -s "${dir}/scripts/mosca.py" "${PREFIX}/bin/"

# this file is not in use until bioconda problem is fixed