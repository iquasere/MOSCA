# -*- coding: utf-8 -*-
"""
MOSCA's Metaproteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

import glob
import os
import pathlib
import shlex
import shutil
import subprocess
import pandas as pd
import argparse
import time
from lxml import etree
from mosca_tools import parse_fasta, run_command, sort_alphanumeric, parse_blast, run_pipe_command
from progressbar import ProgressBar
from Bio import Entrez
from urllib.error import HTTPError
import requests


# TODO - integrate apt-get install -y libpwiz-tools poppler-utils


class MetaproteomicsAnalyser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's metaproteomics analysis")
        parser.add_argument("-sf", "--spectra-folder", type=str, help="Folder with spectra to be analysed")
        parser.add_argument("-fl", "--files-list", type=str, help="File listing all input spectra")
        parser.add_argument("-t", "--threads", type=str, help="Number of threads to use.")
        parser.add_argument("-o", "--output", type=str,
                            help="Project directory. Output goes to {output}/Metaproteomics")
        parser.add_argument("-w", "--workflow", type=str, help="Workflow to use", choices=['maxquant', 'compomics'])
        parser.add_argument("-db", "--database", type=str,
                            help="Database file (FASTA format) from metagenomics for protein identification")
        parser.add_argument("-cdb", "--contaminants-database", type=str, default=None,
                            help="Database file (FASTA format) with contaminant sequences")
        parser.add_argument("-mpar", "--metaphlan-result", type=str, help="Results from MetaPhlan taxonomic annotation")
        parser.add_argument("-rtl", "--references-taxa-level", type=str, default='genus',
                            choices=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
                            help="Taxonomic level to retrieve reference proteomes from")
        parser.add_argument("-bs", "--batch-size", type=int, default=5000,
                            help="How many IDs to submit per request to NCBI")
        parser.add_argument("-ma", "--max-attempts", type=str, default=3, help="Maximum attempts to access NCBI")
        parser.add_argument("-exps", "--experiment-names", type=str, help="Names of experiences (comma-separated)")
        parser.add_argument("--protease", type=str, help="Filename in fasta format of protease sequence",
                            default='Trypsin')
        parser.add_argument("-e", "--experiments", type=str, help="Experiments file")
        parser.add_argument("-if", "--input-format", type=str, help="If experiments is in either TSV or excel format",
                            choices=['tsv', 'excel'])
        parser.add_argument("-mmem", "--max-memory", type=str, default=4,
                            help="Maximum memory to use for Peptide-to-Spectrum matching (in Gb)")
        parser.add_argument("-rd", "--resources-directory", type=str, default=os.path.expanduser('~/resources'),
                            help="Directory for storing databases and other important files")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        args.memory *= 1024
        return args

    def get_proteome_uniprot(self, taxid, output):
        res = requests.get(f'https://www.uniprot.org/uniprot/?query=taxonomy:{taxid}&format=fasta')
        with open(output, 'w') as f:
            f.write(res.content.decode('utf8'))

    '''
    Input:
        mpa_result: str - filename of MetaPhlan result
        output: str - name of folder to output reference proteomes
        references_taxa_level: str - one of 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
        batch_size: int - number of IDs per request to NCBI
        max_attemps: int - number of requests failed before giving up
    '''

    def add_reference_proteomes(self, mpa_result, output, references_taxa_level='genus'):
        taxa_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        mpa_data = pd.read_csv(mpa_result, sep='\t', skiprows=3)
        mpa_data = mpa_data[mpa_data['NCBI_tax_id'].str.count('\|') == taxa_levels.index(references_taxa_level)]
        taxids = [ide.split('|')[-1] for ide in mpa_data['NCBI_tax_id']]
        print(f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt.')
        pbar = ProgressBar()
        for taxid in pbar(taxids):
            self.get_proteome_uniprot(taxid, '{}/{}.fasta'.format(output, taxid))
        return taxids

    '''   
    input: 
        database: str - name of fasta file from metagenomics
        output: str - output folder for protease file and final database
        contaminants_database: str - filename of cRAP database
        protease: str - protease used - still only trypsin
    output:
        a FASTA file named output will be created, containing the sequences from reference proteomes, metagenomics, 
        contaminants database and of the protease
    '''

    def database_generation(self, database, output, mpa_result, contaminants_database=None, protease='trypsin',
                            references_taxa_level='genus', batch_size=5000, max_attempts=3, threads='1'):
        print(f'Generating new database in {output}')

        # Get reference proteomes for the various taxa
        taxids = self.add_reference_proteomes(mpa_result, output, references_taxa_level,
                                              batch_size=batch_size, max_attempts=max_attempts)

        # Add protease
        if protease == 'Trypsin':
            if not os.path.isfile('{}/P00761.fasta'.format(output)):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                run_command(f'wget https://www.uniprot.org/uniprot/P00761.fasta -P {output}')
            protease = 'f{output}/P00761.fasta'
        else:  # is an inputed file
            if not os.path.isfile(protease):
                exit(f'Protease file does not exist: {protease}')

        files = [f'{output}/{taxid}.fasta' for taxid in taxids] + [database, protease]

        if contaminants_database is not None:
            self.verify_crap_db(contaminants_database)
            files.append(contaminants_database)

        run_command(f"cat {' '.join(files)}", output=f'{output}/predatabase.fasta', mode='w')

        # Join aminoacid lines, and remove empty lines
        run_pipe_command(f"awk '{{if ($0 ~ /^>/) {{print \"\\n\"$0}} else {{printf $0}}}}' {output}/predatabase.fasta",
                         output=f'{output}/database.fasta')

        run_command(f'seqkit rmdup -s -i -w 0 -o {output}/unique.fasta -D {output}/seqkit_duplicated.detail.txt '
                    f'-j {threads} {output}/database.fasta')

        # Remove asterisks (non identified aminoacids) and plicas
        run_pipe_command(f"""sed -i "s/[*\']//g" {output}/unique.fasta""")

        run_command(f'seqkit rmdup -n -i -w 0 -o {output}/1st_search_database.fasta '
                    f'-D {output}/seqkit_duplicated.detail.txt -j {threads} {output}/unique.fasta')

    '''
    Input:
        files_list: str - filename of list of input spectra files
        out_dir: str - foldername of output
        format: str - format to convert to the spectra (mgf,)
    '''

    def convert_spectra_format(self, files_list, out_dir, format='mgf', peak_picking=False):
        run_pipe_command(f"""msconvert -f {files_list} --{format} -o {out_dir}
                        {' --filter "peakPicking cwt"' if peak_picking else ''}""")

    '''   
    input: 
        crap_folder: folder where crap.fasta should be
    output: 
        confirmation of cRAP database presence in the folder, or download of the 
        database if absent
    '''

    def verify_crap_db(self, contaminants_database='MOSCA/Databases/metaproteomics/crap.fasta'):
        if os.path.isfile(contaminants_database):
            print(f'cRAP database exists at {contaminants_database}')
        else:
            print(f'cRAP database not found at {contaminants_database}. Downloading cRAP database.')
            run_command(f'wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta -O {contaminants_database}')

    '''   
    input:
        protein_fasta: fasta file with proteins from MG analysis, plus trypsin 
        and cRAP sequences
    output:
        a FASTA file named "database + _concatenated_target_decoy.fasta" will be 
        created, containing interleaved original and decoy sequences
    '''

    def create_decoy_database(self, database):
        decoy_database = database.replace('.fasta', '_concatenated_target_decoy.fasta')
        if not os.path.isfile(decoy_database):
            run_command(f'searchgui eu.isas.searchgui.cmd.FastaCLI -in {database} -decoy')
        else:
            print(f'{decoy_database} already exists!')

    '''   
    input: 
        output: name of parameters file
        database: name of FASTA decoy database
        protein_fdr: float - FDR at the protein level in percent
    output:
        a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
    '''

    def generate_parameters_file(self, output, database, protein_fdr=1):
        run_pipe_command(
            f'searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -out {output} -db {database} -prec_tol 10 '
            f'-frag_tol 0.02 -enzyme Trypsin -fixed_mods "Carbamidomethylation of C" -variable_mods '
            f'"Oxidation of M, Acetylation of protein N-term" -mc 2 -protein_fdr {protein_fdr}')

    '''   
    input: 
        spectra_folder: folder containing the raw spectra files
        output: folder to output results
        parameters_file: parameters filename
        search_engines: search engines to perform PSM
    output:
        a "searchgui_out.zip" file will be created in the output folder
    '''

    def peptide_spectrum_matching(self, spectra_folder, output, parameters_file, threads='12',
                                  search_engines=('xtandem', 'myrimatch', 'msgf'), max_memory=4096):
        run_command(
            f"searchgui eu.isas.searchgui.cmd.SearchCLI -Xmx{max_memory}M -spectrum_files {spectra_folder} "
            f"-output_folder {output} -id_params {parameters_file} -threads {threads}"
            f"{''.join([f' -{engine} 1' for engine in search_engines])}")

    '''   
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
    '''

    def browse_identification_results(self, spectra_folder, parameters_file,
                                      searchcli_output, peptideshaker_output, experiment_name='experiment',
                                      sample_name='sample', replicate_number='1', max_memory=4096):
        try:
            run_command(
                f'peptide-shaker -Xmx{max_memory}M eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files '
                f'{spectra_folder} -experiment {experiment_name} -sample {sample_name} -replicate {replicate_number} '
                f'-identification_files {searchcli_output} -out {peptideshaker_output}')
        except:
            print('Producing Peptide-Shaker result failed! Maybe no identifications were obtained?')

    '''   
    input: 
        peptideshaker_output: peptideshaker output filename
        reports_folder: folder to where output reports
        reports_list: list of INTEGERS from 0 to 11 corresponding to the reports
        to output
    output:
        if it doesn't exist, "reports_folder" will be created
        reports will be outputed to "reports_folder"
    '''

    def generate_reports(self, peptideshaker_output, reports_folder, reports_list=[str(n) for n in range(12)]):
        print(f'Created {reports_folder}')
        pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)  # creates folder for reports
        run_command(
            f"peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in {peptideshaker_output} "
            f"-out_reports {reports_folder} -reports {','.join(reports_list)}")

    '''
    Input:
    Output:
    '''

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

    '''   
    input: 
        protein_report: name of file containing protein report from PeptideShaker
        output: name of file to output
        blast: name of blast file with UniProt annotations - only needed if uniprot_ids = False
        uniprot_ids: boolean, True if protein sequences already are named with UniProt IDs
            If False, will need blast file to retrieve the UniProt IDs from
        samples_names: names of samples respective to protein_reports
    output:
        A tab separated spectra count file named "output", with UniProt IDs
        and corresponding spectra quantification
    '''

    def spectra_counting(self, protein_reports, output, blast=None, uniprot_ids=False,
                         samples_names=None):
        protein_reports = sort_alphanumeric(protein_reports)
        if samples_names is None:
            samples_names = [filename.split('/')[-3] for filename in
                             protein_reports]  # samples names are the folder containing the reports folder
        spectra_count = pd.DataFrame(columns=['Main Accession'])
        for i in range(len(protein_reports)):
            report = pd.read_csv(protein_reports[i], sep='\t', index_col=0)
            if not uniprot_ids:
                blast = parse_blast(blast)
                blast['sseqid'] = [ide.split('|')[-1] for ide in blast.sseqid]
                report = pd.merge(report, blast[['qseqid', 'sseqid']], left_on='Main Accession', right_on='qseqid')
                report = report[['sseqid', '#PSMs']]
            else:
                report = report[['Main Accession', '#PSMs']]
            report.columns = ['Main Accession', samples_names[i]]
            report = report.groupby('Main Accession')[samples_names[i]].sum().reset_index()
            spectra_count = pd.merge(spectra_count, report, on='Main Accession', how='outer')
        spectra_count[samples_names] = spectra_count[samples_names].fillna(value=0).astype(int)
        spectra_count.to_csv(output, sep='\t', index=False)

    '''
    Metaproteomics with MaxQuant
    '''

    '''
    input:
        output: name of standard parameters file to create
    output:
        a standard parameters file for MaxQuant named "output" will be created
    '''

    def create_mqpar(self, output):
        if os.path.isfile(output):  # the create file command will not create a new one if the file already exists
            os.remove(
                output)  # even if that file already has non-default information in it, messing with the next commands
        run_command(f'maxquant {output} --create')

    '''
    input:
        mqpar: name of the mqpar.xml file to have its values changed
        fasta_database: name of the FASTA database for metaproteomics
        spectra_folder: name of the folder with RAW spectra files
        experiment_names: list of experiments, with one element for each file
    output:
        the "file" file will be updated with the new parameters
    '''

    def edit_maxquant_mqpar(self, mqpar, fasta_database, spectra_folder, experiment_names, threads=1,
                            spectra_format='RAW', protein_fdr=0.01):
        print('Updating parameters file information.')
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(mqpar, parser)
        root = tree.getroot()
        root.find("fastaFiles/FastaFileInfo/fastaFilePath").text = fasta_database
        print('Fasta database = ' + fasta_database)
        root.find("separateLfq").text = 'True'
        root.find("numThreads").text = str(threads)
        root.find("proteinFdr").text = str(protein_fdr)
        print('Number of threads = ' + str(threads))
        for child in ['filePaths/string', 'experiments/string', 'fractions/short', 'ptms/boolean',
                      'paramGroupIndices/int']:
            tree.xpath(child)[0].getparent().remove(tree.xpath(child)[
                                                        0])  # the default params file of MaxQuant brings default arguments that can't be present
        filePaths = root.find("filePaths")
        experiments = root.find("experiments")
        fractions = root.find("fractions")
        ptms = root.find("ptms")
        paramGroupIndices = root.find("paramGroupIndices")
        files = sort_alphanumeric(glob.glob(f'{spectra_folder}/*.{spectra_format}'))
        for i in range(len(files)):
            print(f'Adding file: {files[i]}')
            etree.SubElement(filePaths, 'string').text = files[i]
            print(f'Experiment = {experiment_names[i]}')
            etree.SubElement(experiments, 'string').text = experiment_names[i]
            etree.SubElement(fractions, 'short').text = '32767'
            etree.SubElement(ptms, 'boolean').text = 'True'
            etree.SubElement(paramGroupIndices, 'int').text = '0'
        root.find("parameterGroups/parameterGroup/lfqMode").text = '1'
        tree.write(mqpar, pretty_print=True)
        print(f'Parameters file is available at {mqpar}')

    '''
    input:
        mqpar: name of MaxQuant parameters file
        spectra_folder: folder containing the spectra RAW files
        output_folder: folder where the "combined" folder will be placed
    output:
        
    '''

    def run_maxquant(self, mqpar, spectra_folder, output_folder):
        for directory in [f'{spectra_folder}/combined', output_folder]:
            if os.path.isdir(directory):
                shutil.rmtree(directory, ignore_errors=True)
        run_command(f'maxquant {mqpar}')  # TODO - get the shell messages from MaxQuant to appear
        os.rename(f'{spectra_folder}/combined', output_folder)

    '''
    Input:
        database: str - filename of database for protein identification
        output: str - name of folder to store results
        spectra_folder: str - name of folder with input spectra
        threads: int - number of threads to use
        sample_name: str
        experiment_name: str
        replicate_number: str
        protein_fdr: int - FDR at the protein level in percent
    '''

    def compomics_workflow(self, database, output, spectra_folder, threads=1, sample_name='sample',
                           experiment_name='experiment', replicate_number='1', protein_fdr=1, max_memory=4096):

        self.create_decoy_database(database)

        try:  # try/except - https://github.com/compomics/searchgui/issues/217
            self.generate_parameters_file(
                f'{output}/params.par',
                database.replace('.fasta', '_concatenated_target_decoy.fasta'),
                protein_fdr=protein_fdr)
        except:
            print('An illegal reflective access operation has occurred. But MOSCA can handle it.')

        self.peptide_spectrum_matching(spectra_folder, output, f'{output}/params.par', threads=threads,
                                       max_memory=max_memory)

        self.browse_identification_results(
            spectra_folder,
            f'{output}/params.par',
            f'{output}/searchgui_out.zip',
            f'{output}/ps_output.cpsx',
            max_memory=max_memory,
            sample_name=sample_name,
            experiment_name=experiment_name,
            replicate_number=replicate_number)

        try:  # try/except - if no identifications are present, will throw an error
            self.generate_reports(
                f'{output}/ps_output.cpsx',
                f'{output}/reports')
        except:
            print('No identifications?')

        # self.spectra_counting('{}/reports/{}_{}_{}_Default_Protein_Report.txt'.format(
        #    output, experiment_name, sample_name, replicate_number), self.blast, self.output + '/Spectra_counting.tsv')

    '''
    Input:
        mqpar: str - filename of parameters file to create
        database: str - filename of protein database
        spectra_folder: str - name of folder containing spectra
        experiment_names: list - of names to use as experiment names
        output: str - name of folder to output maxquant results
        threads: int - number of threads to use
    '''

    def maxquant_workflow(self, mqpar, database, spectra_folder, experiment_names, output, threads=1,
                          spectra_format='RAW', protein_fdr=0.01):
        self.create_mqpar(mqpar)
        self.edit_maxquant_mqpar(mqpar, database, spectra_folder, experiment_names, threads=threads,
                                 spectra_format=spectra_format, protein_fdr=protein_fdr)
        self.run_maxquant(mqpar, spectra_folder, output)

    def select_proteins_for_second_search(self, original_db, output, results_files, column='Main Accession'):
        proteins = list()
        for file in results_files:
            proteins += pd.read_csv(file, sep='\t')[column].tolist()
        proteins = set(proteins)
        print(f'Selected {len(proteins)} proteins for 2nd peptide-to-spectrum matching.')
        with open(f'{output}/2nd_search_ids.txt', 'w') as f:
            f.write('\n'.join(proteins))
        run_pipe_command(
            f"seqkit grep {original_db} -w 0 --id-regexp '\|(.*)\|' -f {output}/2nd_search_ids.txt",
            output=f'{output}/2nd_search_database.fasta'
        )  # the '\|(.*)\|' rule also works for the selection of full ids (like the ones coming from FragGeneScan)

    def run(self):

        args = self.get_arguments()
        experiments = (pd.read_csv(args.experiments, sep='\t') if args.input_format == 'tsv' else
                       pd.read_excel(args.experiments))
        experiments = experiments[experiments['Data type'] == 'protein']

        for sample in set(experiments['Sample']):
            partial_sample = experiments[experiments['Sample'] == sample].reset_index(drop=True)

            for foldername in ['spectra', '1st_search', '2nd_search']:
                pathlib.Path(f'{args.output}/Metaproteomics/{sample}/{foldername}').mkdir(parents=True, exist_ok=True)

            self.database_generation(
                args.database, f'{args.output}/Metaproteomics/{sample}', args.metaphlan_result,
                contaminants_database=args.contaminants_database, protease=args.protease,
                references_taxa_level=args.references_taxa_level, batch_size=args.batch_size,
                max_attempts=args.max_attemps)

            if args.workflow == 'maxquant':
                self.maxquant_workflow(
                    f'{args.output}/mqpar.xml', f'{args.output}/database.fasta', args.spectra_folder,
                    args.experiment_names.split(','), args.output, threads=args.threads, spectra_format='RAW',
                    protein_fdr=1)

                # TODO - have to fix this next command - put right names and parameters
                self.maxquant_workflow('{}/mqpar.xml'.format(args.output), '{}/database.fasta'.format(args.output),
                                       args.spectra_folder, args.experiment_names.split(','), args.output,
                                       threads=args.threads, spectra_format='RAW', protein_fdr=0.01)

            elif args.workflow == 'compomics':
                for name in set(partial_sample['Name']):
                    if os.path.isfile(
                            '{0}/Metaproteomics/{1}/1st_search/{2}/reports/{1}_{2}_1_Default_Protein_Report_with_non-validated_matches.txt'.format(
                                args.output, sample, name)):
                        print(f'[{name}] already analysed!')
                        continue

                    for foldername in ['spectra', '1st_search', '2nd_search']:
                        pathlib.Path('{}/Metaproteomics/{}/{}/{}'.format(args.output, sample, foldername, name)).mkdir(
                            parents=True, exist_ok=True)
                    partial_name = partial_sample[partial_sample['Name'] == name].reset_index(drop=True)
                    files2convert = list()

                    for i in range(len(partial_name)):
                        if partial_name.iloc[i]['Files'].endswith('_0.mzML'):
                            continue
                        if not partial_name.iloc[i]['Files'].endswith('.mgf'):
                            files2convert.append(partial_name.iloc[i]['Files'])
                        else:
                            shutil.copy(partial_name.iloc[i]['Files'], '{}/Metaproteomics/{}/spectra/{}'.format(
                                args.output, sample, name))
                    # Conversion to Mascot Generic Format
                    if len(files2convert) > 0:
                        with open('{}/Metaproteomics/{}/files2convert.txt'.format(args.output, sample), 'w') as f:
                            f.write('\n'.join(files2convert))
                        self.convert_spectra_format('{}/Metaproteomics/{}/files2convert.txt'.format(args.output, sample),
                                                    '{}/Metaproteomics/{}/spectra/{}'.format(args.output, sample, name),
                                                    format='mgf')
                    
                    self.compomics_workflow('{}/Metaproteomics/{}/1st_search_database.fasta'.format(args.output, sample),
                                            '{}/Metaproteomics/{}/1st_search/{}'.format(args.output, sample, name),
                                            '{}/Metaproteomics/{}/spectra/{}'.format(args.output, sample, name),
                                            threads=args.threads, protein_fdr=100, max_memory=args.max_memory,
                                            experiment_name=sample, sample_name=name)
                    
                self.select_proteins_for_second_search(
                    '{}/Metaproteomics/{}/1st_search_database.fasta'.format(args.output, sample),
                    '{}/Metaproteomics/{}'.format(args.output, sample),
                    glob.glob('{}/Metaproteomics/{}/1st_search/*/reports/*_Protein_Report_with_non-validated_matches.txt'.format(
                        args.output, sample)), column='Main Accession')

                for name in set(partial_sample['Name']):
                    self.compomics_workflow(
                        '{}/Metaproteomics/{}/2nd_search_database.fasta'.format(args.output, sample),
                        '{}/Metaproteomics/{}/2nd_search/{}'.format(args.output, sample, name),
                        '{}/Metaproteomics/{}/spectra/{}'.format(args.output, sample, name),
                        threads=args.threads, protein_fdr=1, max_memory=args.max_memory,
                        experiment_name=sample, sample_name=name)

                self.join_ps_reports(
                    ['{0}/Metaproteomics/{1}/2nd_search/{2}/reports/{1}_{2}_1_Default_Protein_Report.txt'.format(
                        args.output, sample, name) for name in set(partial_sample['Name'])],
                    local_fdr=5, validation=False
                ).to_csv('{}/Metaproteomics/{}/quantification_table.tsv'.format(args.output, sample), sep='\t')

            else:
                print('Not a valid workflow option!')
            exit()

    def background_inputation(self, df):
        return df.fillna(value=df.min().min())

    def censored_inputation(self, df, replicates):
        for replicate in replicates:
            for line in df.index:
                if df.loc[df.index[0]][replicates[0]].isnull().sum() > 1:
                    df.loc[df.index[0]][replicates[0]] = (
                        self.background_inputation(df.loc[df.index[0]][replicates[0]]))


def censored_inputation(df, replicates):
    pbar = ProgressBar()
    min_val = df.min().min()
    for replicate in pbar(replicates):
        for line in df.index:
            if df.loc[line][replicate].isnull().sum() > 1:
                df.loc[line][replicate] = df.loc[line][replicate].fillna(value=min_val)
    return df


def background_inputation(df):
    return df.fillna(value=df.min())


if __name__ == '__main__':
    MetaproteomicsAnalyser().run()
