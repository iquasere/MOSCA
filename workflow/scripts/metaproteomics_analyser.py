# -*- coding: utf-8 -*-
"""
MOSCA's Metaproteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

from glob import glob
import os
import pathlib
import shutil
import pandas as pd
import argparse
from tqdm import tqdm
import requests
from time import sleep
from mosca_tools import run_command, run_pipe_command, multiprocess_fun


class MetaproteomicsAnalyser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's metaproteomics analysis")
        parser.add_argument("-sfs", "--spectra-folders", help="List of folders with spectra to be analysed")
        parser.add_argument("-ns", "--names", help="List of names for datasets")
        parser.add_argument("-t", "--threads", help="Number of threads to use.")
        parser.add_argument("-o", "--output", help="Output directory, where results are stored.")
        parser.add_argument(
            "-db", "--database", help="Database file (FASTA format) from metagenomics for protein identification")
        parser.add_argument(
            "-cdb", "--contaminants-database", default=None,
            help="Database file (FASTA format) with contaminant sequences")
        parser.add_argument("-ur", "--upimapi-result", help="Results from UPIMAPI")
        parser.add_argument(
            "-bs", "--batch-size", type=int, default=5000, help="How many IDs to submit per request to NCBI")
        parser.add_argument("-ma", "--max-attempts", default=3, help="Maximum attempts to access NCBI")
        parser.add_argument("--protease", help="Filename in fasta format of protease sequence", default='Trypsin')
        parser.add_argument(
            "-mmem", "--max-memory", type=int, default=4, help="Maximum memory (Gb) to use for Peptide-to-Spectrum matching [4]")
        parser.add_argument(
            "-rd", "--resources-directory", default=os.path.expanduser('~/resources'),
            help="Directory for storing databases and other important files")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        args.spectra_folders = args.spectra_folders.split(',')
        args.names = args.names.split(',')
        if len(args.names) != len(args.spectra_folders):
            exit('Length of names and spectra_folders must be equal.')
        args.max_memory *= 1024
        return args

    def get_proteome_uniprot(self, taxid, max_tries=3):
        tries = 0
        done = False
        while tries < max_tries and not done:
            try:
                res = requests.get(
                    f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28taxonomy_id%3A{taxid}%29')
                done = True
            except:
                print(f'Failed! {max_tries - tries} tries remaining.')
                tries += 1
                sleep(10)
        return res.content.decode('utf8')

    def add_reference_proteomes(self, upimapi_res, output):
        taxids = pd.read_csv(
            upimapi_res, sep='\t')['Taxonomic lineage IDs (SPECIES)'].dropna().unique().astype(int).tolist()
        with open(output, 'w') as f:
            for taxid in tqdm(taxids, desc=f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt'):
                f.write(self.get_proteome_uniprot(taxid))

    def database_generation(
            self, mg_orfs, output, upimapi_res, contaminants_database=None, protease='Trypsin', threads=1):
        """
        Build database from MG analysis
        :param mg_orfs:
        :param output:
        :param mpa_result:
        :param contaminants_database:
        :param protease:
        :param threads:
        """
        print(f'Generating new database in {output}')
        # Get reference proteomes for the various taxa
        self.add_reference_proteomes(upimapi_res, f'{output}/ref_proteomes.fasta')
        # Add protease
        if protease == 'Trypsin':
            if not os.path.isfile(f'{output}/P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                run_command(f'wget https://www.uniprot.org/uniprot/P00761.fasta -P {output}')
            protease = f'{output}/P00761.fasta'
        else:  # is an inputed file
            if not os.path.isfile(protease):
                exit(f'Protease file does not exist: {protease}')
        files = [mg_orfs, f'{output}/ref_proteomes.fasta', protease]
        if contaminants_database is not None:
            self.verify_crap_db(contaminants_database)
            files.append(contaminants_database)
        run_command(f"cat {' '.join(files)}", output=f'{output}/predatabase.fasta', mode='w')
        # Join aminoacid lines, and remove empty lines
        run_pipe_command(f"awk '{{if ($0 ~ /^>/) {{print \"\\n\"$0}} else {{printf $0}}}}' {output}/predatabase.fasta",
                         output=f'{output}/database.fasta')
        run_pipe_command(f"awk '{{print $1}}' {output}/database.fasta", output=f'{output}/predatabase.fasta')
        # Remove asterisks (non identified aminoacids) and plicas
        run_pipe_command(f"""sed -i "s/[*\']//g" {output}/predatabase.fasta""")
        # Remove duplicate sequences
        run_command(
            f'seqkit rmdup -s -i -w 0 -o {output}/1st_search_database.fasta -D {output}/seqkit_duplicated.detail.txt '
            f'-j {threads} {output}/predatabase.fasta')
        for file in ['database.fasta', 'predatabase.fasta', 'ref_proteomes.fasta']:
            os.remove(f'{output}/{file}')

    def raw_to_mgf(self, file, out_dir):
        """
        Convert raw file to mgf
        :param file:
        :param out_dir:
        :param outfmt:
        :param peak_picking:
        :return:
        """
        folder, filename = os.path.split(file)
        pathlib.Path('named_volume').mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f'{folder}/{filename}', f'/data/{filename}')
        run_pipe_command(
            f'docker run --rm -e WINEDEBUG=-all -v named_volume:/data '
            f'chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/{filename} --mgf '
            f'--filter "peakPicking cwt"')
        shutil.copyfile(
            f'/data/{".".join(filename.split(".")[:-1])}.mgf', f'{out_dir}/{".".join(filename.split(".")[:-1])}.mgf')

    def spectra_in_proper_state(self, folder, out_dir, threads=1):
        """
        Convert raw spectra files to MGF, and put all spectra in out_dir
        :param threads:
        :param folder:
        :param out_dir:
        :return:
        """
        already_mgfs, files2convert = [], []
        for file in glob(f'{folder}/*'):
            if os.path.isfile(file):
                if file.endswith('.mgf'):
                    already_mgfs.append(file)
                else:
                    files2convert.append(file)
        if len(files2convert) > 0:
            multiprocess_fun(self.raw_to_mgf, [(file, out_dir) for file in files2convert], threads=threads)
        for file in already_mgfs:
            shutil.copyfile(f'{folder}/{file}', f'{out_dir}/{file}')

    def verify_crap_db(self, contaminants_database='MOSCA/Databases/metaproteomics/crap.fasta'):
        """
        Checks if contaminants database exists. If not, cRAP will be downloaded.
        :param contaminants_database:
        :return:
        """
        if os.path.isfile(contaminants_database):
            print(f'cRAP database exists at {contaminants_database}')
        else:
            print(f'cRAP database not found at {contaminants_database}. Downloading cRAP database.')
            run_command(f'wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta -O {contaminants_database}')

    def create_decoy_database(self, database):
        """
        input:
            protein_fasta: fasta file with proteins from MG analysis, plus trypsin
            and cRAP sequences
        output:
            a FASTA file named "database + _concatenated_target_decoy.fasta" will be
            created, containing interleaved original and decoy sequences
        """
        decoy_database = database.replace('.fasta', '_concatenated_target_decoy.fasta')
        if not os.path.isfile(decoy_database):
            run_command(f'searchgui eu.isas.searchgui.cmd.FastaCLI -in {database} -decoy')
        else:
            print(f'{decoy_database} already exists!')

    def generate_parameters_file(self, output, protein_fdr=1):
        """
        input:
            output: name of parameters file
            database: name of FASTA decoy database
            protein_fdr: float - FDR at the protein level in percent
        output:
            a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
        """
        run_pipe_command(
            f'searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} -prec_tol 10 '
            f'-frag_tol 0.02 -enzyme Trypsin -fixed_mods "Carbamidomethylation of C" -variable_mods '
            f'"Oxidation of M, Acetylation of protein N-term" -mc 2 -protein_fdr {protein_fdr}')

    def split_database(self, database, n_proteins=5000000):
        """
        Split database into smaller files
        :param database:
        :param n_proteins:
        :return:
        """
        run_command(f'seqkit split -s {n_proteins} {database} -O {os.path.dirname(database)}')

    def peptide_spectrum_matching(
            self, spectra_folder, output, parameters_file, database, threads=12, max_memory=4096,
            search_engines=('xtandem', 'myrimatch', 'msgf')):
        """
        input:
            spectra_folder: folder containing the raw spectra files
            output: folder to output results
            parameters_file: parameters filename
            search_engines: search engines to perform PSM
        output:
            a "searchgui_out.zip" file will be created in the output folder
        """
        run_command(
            f"searchgui eu.isas.searchgui.cmd.SearchCLI -Xmx{max_memory}M -spectrum_files {spectra_folder} "
            f"-output_folder {output} -id_params {parameters_file} -threads {threads} -fasta_file {database} "
            f"{' '.join([f'-{engine} 1' for engine in search_engines])}")

    def browse_identification_results(
            self, spectra_folder, output, parameters_file, database, searchcli_output,  name, max_memory=4096):
        """
        input:
            spectra_folder: folder containing the raw spectra files
            output: folder to output results
            parameters_file: parameters filename
            searchcli_output: searchcli output filename
            peptideshaker_output: peptideshaker output filename
            experiment_name: name of experiment
            sample_name: name of sample
            replicate_number: number of replicate (STRING)
        output:
            a file will be outputed with the validation of the PSMs and protein
            identifications
        """
        try:
            run_command(
                f'peptide-shaker -Xmx{max_memory}M eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files '
                f'{spectra_folder} -reference {name} -identification_files {searchcli_output} '
                f'-out {output} -fasta_file {database} -id_params {parameters_file}')
        except:
            print('Producing Peptide-Shaker result failed! Maybe no identifications were obtained?')

    def generate_reports(self, peptideshaker_output, reports_folder, reports=["10"], max_memory=40960):
        """
        input:
            peptideshaker_output: peptideshaker output filename
            reports_folder: folder to where output reports
            reports_list: list of INTEGERS from 0 to 11 corresponding to the reports
            to output
        output:
            if it doesn't exist, "reports_folder" will be created
            reports will be outputed to "reports_folder"
        """
        print(f'Created {reports_folder}')
        pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)  # creates folder for reports
        run_command(
            f"peptide-shaker -Xmx{max_memory}M eu.isas.peptideshaker.cmd.ReportCLI -in {peptideshaker_output} "
            f"-out_reports {reports_folder} -reports {','.join(reports)}")

    def join_ps_reports(self, files, local_fdr=None, validation=False):
        result = pd.DataFrame(columns=['Main Accession'])
        for file in files:
            data = pd.read_csv(file, sep='\t', index_col=0)
            if local_fdr is not None:
                data = data[data['Confidence [%]'] > 100 - local_fdr]
            if validation:
                data = data[data['Validation'] == 'Confident']
            if len(data) > 0:
                name = file.split('/')[-1].split('_')[1]
                data = data[['Main Accession', '#PSMs']]
                data.columns = ['Main Accession', name]
                result = pd.merge(result, data, on='Main Accession', how='outer')
        return result

    def compomics_run(
            self, database, output, spectra_folders, name, params, threads=1, max_memory=4096, reports=['10']):
        """
        Run compomics workflow on the given spectra folders
        :param params:
        :param database:
        :param output:
        :param spectra_folders:
        :param name:
        :param threads:
        :param max_memory:
        :param reports:
        :return:
        """
        self.peptide_spectrum_matching(
            spectra_folders, output, params, database, threads=threads, max_memory=max_memory)
        self.browse_identification_results(
            spectra_folders, f'{output}/ps_output.psdb', params, database, f'{output}/searchgui_out.zip',
            name, max_memory=max_memory)
        try:
            self.generate_reports(
                f'{output}/ps_output.psdb', f'{output}/reports', max_memory=max_memory, reports=reports)
        except Exception as e:
            print(f'Reporter generation sent an error: {e}\nBut likely worked!')

    def select_proteins_for_second_search(self, original_db, output, results_files, column='Main Accession'):
        proteins = []
        for file in results_files:
            proteins += pd.read_csv(file, sep='\t')[column].tolist()
        proteins = set(proteins)
        print(f'Selected {len(proteins)} proteins for 2nd peptide-to-spectrum matching.')
        with open(f'{output}/2nd_search_ids.txt', 'w') as f:
            f.write('\n'.join(proteins))
        # the '\|(.*)\|' rule also works for the selection of full ids (like the ones coming from FragGeneScan)
        run_pipe_command(
            f"seqkit grep {original_db} -w 0 --id-regexp '\|(.*)\|' -f {output}/2nd_search_ids.txt",
            output=f'{output}/2nd_search_database.fasta')

    def run(self):
        args = self.get_arguments()
        '''
        self.database_generation(
            args.database, args.output, args.upimapi_result, contaminants_database=args.contaminants_database,
            protease=args.protease)
        self.create_decoy_database(f'{args.output}/1st_search_database.fasta')
        self.split_database(
            f'{args.output}/1st_search_database_concatenated_target_decoy.fasta', n_proteins=5000000)
        try:  # try/except - https://github.com/compomics/searchgui/issues/217
            self.generate_parameters_file(f'{args.output}/1st_params.par', protein_fdr=100)
        except:
            print('An illegal reflective access operation has occurred. But MOSCA can handle it.')
        '''
        proteins_for_second_search = []

        for i in range(len(args.names)):
            out = f'{args.output}/{args.names[i]}'
            for foldername in ['spectra', '2nd_search']:
                pathlib.Path(f'{out}/{foldername}').mkdir(parents=True, exist_ok=True)
            self.spectra_in_proper_state(args.spectra_folders[i], f'{out}/spectra')
            i = 0
            for database in sorted(glob(f'{os.path.dirname(args.output)}/1st_search_database_*.part_*.fasta')):
                i += 1
                pathlib.Path(f'{out}_{i}').mkdir(parents=True, exist_ok=True)
                self.compomics_run(
                    database, f'{out}_{i}/1st_search', f'{out}/spectra', args.names[i],
                    f'{args.output}/1s_params.par', threads=args.threads, protein_fdr=100,
                    max_memory=args.max_memory, reports=['4'])
                df = pd.read_csv(
                    f'{out}_{i}/reports/{args.names[i]}_Default_PSM_Report_with_non-validated_matches.txt',
                    sep='\t', index_col=0)
                for protein_group in df['Protein(s)'].str.split(','):
                    proteins_for_second_search += protein_group

        with open(f'{args.output}/2nd_search_proteins.txt', 'w') as f:
            f.write('\n'.join(set(proteins_for_second_search)))
        run_command(
            f'seqkit grep -f {args.output}/2nd_search_proteins.txt -o {args.output}/2nd_search_database.fasta '
            f'{args.output}/1st_search_database.fasta')

        self.create_decoy_database(f'{args.output}/2nd_search_database.fasta')
        self.generate_parameters_file(f'{args.output}/2nd_params.par', protein_fdr=1)

        # TODO - check if will be necessary to split the database
        # self.split_database(
        #    f'{args.output}/2nd_search_database_concatenated_target_decoy.fasta', n_proteins=5000000)

        spectracounts = pd.DataFrame(columns=['Main Accession'])
        for name in args.names:
            out = f'{args.output}/{name}'
            self.compomics_run(
                f'{args.output}/2nd_search_database.fasta', f'{out}/2nd_search', f'{out}/spectra', name,
                f'{args.output}/2nd_params.par', threads=args.threads, max_memory=args.max_memory,
                reports=['9', '10'])
            spectracounts = pd.merge(spectracounts, pd.read_csv(
                f'{out}/2nd_search/reports/{name}_Default_Protein_Report.txt', sep='\t', index_col=0
            )[['Main Accession', '#PSMs']].rename(columns={'#PSMs': name}), how='outer', on='Main Accession')
        spectracounts.groupby('Main Accession')[args.names].sum().reset_index().fillna(value=0.0).to_csv(
            f'{args.output}/spectracounts.tsv', sep='\t', index=False)


if __name__ == '__main__':
    MetaproteomicsAnalyser().run()
