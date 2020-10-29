dir="${PREFIX}/share/MOSCA"
mkdir -p "${dir}/scripts" "${PREFIX}/bin"
cp workflow/scripts/* "${dir}/scripts"
cp workflow/Snakefile workflow/mosca.py "${dir}/scripts"
cp -r resources "${dir}/resources"
chmod +x "${dir}/scripts/mosca.py"
ln -s "${dir}/scripts/mosca.py" "${PREFIX}/bin/"

# Databases download
resources_dir="${PREFIX}/share/MOSCA/resources"
svn export https://github.com/timflutre/trimmomatic/trunk/adapters "${resources_dir}/illumina_adapters"
svn export https://github.com/biocore/sortmerna/trunk/rRNA_databases "${resources_dir}/rRNA_databases"
find "${resources_dir}/rRNA_databases/*" | grep -v ".fasta" | xargs rm -fr

# Proteomics tools installation
apt-get install -y libpwiz-tools poppler-utils
perl ~/anaconda3/opt/krona/install.pl
wget http://genesis.ugent.be/maven2/eu/isas/searchgui/SearchGUI/3.3.16/SearchGUI-3.3.16-mac_and_linux.tar.gz
tar -xzf SearchGUI-3.3.16-mac_and_linux.tar.gz
wget http://genesis.ugent.be/maven2/eu/isas/peptideshaker/PeptideShaker/1.16.41/PeptideShaker-1.16.41.zip
unzip PeptideShaker-1.16.41.zip