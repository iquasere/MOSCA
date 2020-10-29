# -*- coding: utf-8 -*-
"""
MOSCA's MetaProteomics class for performing MetaProteomics Analysis

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
from annotation import Annotater
from lxml import etree
from mosca_tools import parse_fasta, run_command, sort_alphanumeric, parse_blast
from progressbar import ProgressBar


class MetaproteomicsAnalyser:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        self.searchgui_exe = glob.glob('SearchGUI-*.*.*/SearchGUI-*.*.*.jar')[-1]
        self.peptide_shaker_exe = glob.glob('PeptideShaker-*.*.*/PeptideShaker-*.*.*.jar')[-1]

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's metaproteomics analysis")
        parser.add_argument("-sf", "--spectra-folder", type=str, required=True,
                            help="Folder with spectra to be analysed")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-o", "--output", type=str, help="Output directory")
        parser.add_argument("-w", "--workflow", type=str, help="Workflow to use", choices=['maxquant','compomics'])
        parser.add_argument("-o", "--output", type=str, help="Output directory")
        parser.add_argument("-rdb", "--reference-database", type=str,
                            help="Database file (FASTA format) with reference sequences for protein identification")
        parser.add_argument("-cdb", "--contaminants-database", type=str,
                            help="Database file (FASTA format) with contaminant sequences")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args
        
    '''   
    input: 
        faa: name of fasta file from FragGeneScan
        output: output folder for protease file and final database
        crap_database: filename of cRAP database
        protease: protease name used
        header: String, defines how the MG information will be handled. May
        significantly reduce the complexity of data and speed up MP analysis
            raw: requires the 'faa' argument. Original information of annotated 
            proteins is not changed, and database will be generated with only 
            the addition of protease and cRAP sequences
            uniprot_sequences: requires the 'blast' argument. Identifiers will be
            retrieved from the blast file, and the protein sequences corresponding 
            to the IDs will be retrieved from UniProt. May take several days!
            uniprot_ids: requires both the 'faa' and 'blast' arguments.
            Protein sequences will be grouped by UniProt ID assigned, and a
            consensus sequence will be generated for each group
    output: 
        a FASTA file named output will be created, containing the sequences from 
        metagenomics, from cRAP database and of the protease
    '''
    def database_generation(self, database, output, crap_database, protease='trypsin', how='raw', blast=None):
        temp = '/'.join(output.split('/')[:-1]) + '/temporary.fasta'
        print('Generating new database in {}'.format(output))

        print('Removing asterisks (*)')
        run_command("sed 's/*//g' {} > {}".format(database, output))

        if protease == 'trypsin':                                                                   # TODO - pig trypsin is not the only protease used, will have to include more - see from searchgui
            if not os.path.isfile('P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                run_command('wget https://www.uniprot.org/uniprot/P00761.fasta -P {}'.format(
                    '/'.join(output.split('/')[:-1])))
            protease = 'P00761.fasta'

        run_command('cat {}'.format(' '.join([database, crap_database, protease])), file=output, mode = 'a')

    '''   
    input: 
        crap_folder: folder where crap.fasta should be
    output: 
        confirmation of cRAP database presence in the folder, or download of the 
        database if absent
    '''    
    def verify_crap_db(self, crap_database = 'MOSCA/Databases/metaproteomics/crap.fasta'):
        if os.path.isfile(crap_database):
            print('cRAP database exists at ' + crap_database)
        else:
            print('cRAP database not found at ' + crap_database + '. Downloading cRAP database.')
            crap_folder = '/'.join(crap_database.split('/')[:-1])
            run_command('wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta -P ' + crap_folder)
            
    '''   
    input:
        protein_fasta: fasta file with proteins from MG analysis, plus trypsin 
        and cRAP sequences
        searchgui_exe: executable .jar of searchgui
    output:
        a FASTA file named "database + _concatenated_target_decoy.fasta" will be 
        created, containing interleaved original and decoy sequences
    ''' 
    def create_decoy_database(self, database, searchgui_exe):
        decoy_database = database.replace('.fasta', '_concatenated_target_decoy.fasta')
        if not os.path.isfile(decoy_database):
            bashCommand = ('java -cp {} eu.isas.searchgui.cmd.FastaCLI -in {} -decoy'.format(
                    searchgui_exe, database))
            run_command(bashCommand)
        else:
            print(decoy_database + ' already exists!')
        
    '''   
    input: 
        output: name of parameters file
        database: name of FASTA decoy database
        searchgui_exe: executable .jar of searchgui
    output:
        a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
    ''' 
    def generate_parameters_file(self, output, database, searchgui_exe):
        bashCommand = ('java -cp {} eu.isas.searchgui.cmd.IdentificationParametersCLI ' + 
                       '-out {} -db {} -prec_tol 10 -frag_tol 0.02 -enzyme Trypsin ' +
                       '-fixed_mods "Carbamidomethylation of C" -variable_mods ' + 
                       '"Oxidation of M, Acetylation of protein N-term" -mc 2').format(
                               searchgui_exe, output, database)
        print(bashCommand)
        bashCommand_correct = shlex.split(bashCommand)                              #the usual way with MoscaTools is not working here, dunno why
        process = subprocess.Popen(bashCommand_correct, stdout=subprocess.PIPE)
        output, error = process.communicate()
        run_command(bashCommand)
        
    '''   
    input: 
        spectra_folder: folder containing the raw spectra files
        output: folder to output results
        parameters_file: parameters filename
        searchgui_exe: executable .jar of searchgui
        search_engines: search engines to perform PSM
    output:
        a "searchgui_out.zip" file will be created in the output folder
    ''' 
    def peptide_spectrum_matching(self, spectra_folder, output, 
                                  parameters_file, searchgui_exe,
                                  search_engines = ['xtandem', 'myrimatch', 'msgf'],
                                  threads = '12'):
        bashCommand = ('java -cp {} eu.isas.searchgui.cmd.SearchCLI -spectrum_files ' +
                       '{} -output_folder {} -id_params {} -threads {}' + 
                       ''.join([' -' + engine + ' 1' for engine in search_engines])).format(
                               searchgui_exe, spectra_folder, output, parameters_file, 
                               threads)
        run_command(bashCommand)
        
    '''   
    input: 
        spectra_folder: folder containing the raw spectra files
        output: folder to output results
        parameters_file: parameters filename
        searchcli_output: searchcli output filename
        peptideshaker_output: peptideshaker output filename
        peptide_shaker_exe: executable .jar of peptide-shaker
        experiment_name: name of experiment
        sample_name: name of sample
        replicate_number: number of replicate (STRING)
    output:
        a file will be outputed with the validation of the PSMs and protein
        identifications
    ''' 
    def browse_identification_results(self, spectra_folder, parameters_file,
                                      searchcli_output, peptideshaker_output,
                                      peptide_shaker_exe, experiment_name = 'MyExperiment', 
                                      sample_name = 'MySample', replicate_number = '1'):
        bashCommand = ('java -cp {} eu.isas.peptideshaker.cmd.PeptideShakerCLI ' +
                       '-spectrum_files {} -experiment {} -sample {} -replicate ' + 
                       '{} -identification_files {} -out {}').format(
                        peptide_shaker_exe, spectra_folder, experiment_name, sample_name, 
                        replicate_number, searchcli_output, peptideshaker_output)
        run_command(bashCommand)
       
    '''   
    input: 
        peptideshaker_output: peptideshaker output filename
        reports_folder: folder to where output reports
        peptide_shaker_exe: executable .jar of peptide-shaker
        reports_list: list of INTEGERS from 0 to 11 corresponding to the reports
        to output
    output:
        if it doesn't exist, "reports_folder" will be created
        reports will be outputed to "reports_folder"
    ''' 
    def generate_reports(self, peptideshaker_output, reports_folder, peptide_shaker_exe,
                         reports_list = [str(n) for n in range(12)]):
        print('Created ' + reports_folder)
        pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)                 # creates folder for reports
        bashCommand = ('java -cp {} eu.isas.peptideshaker.cmd.ReportCLI -in ' + 
                       '{} -out_reports {} -reports ' + ','.join(reports_list)).format(
                               peptide_shaker_exe, peptideshaker_output, reports_folder)
        run_command(bashCommand)
        
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
    def spectra_counting(self, protein_reports, output, blast = None, uniprot_ids = False,
                         samples_names = None):
        protein_reports = sort_alphanumeric(protein_reports)
        if samples_names is None:
            samples_names = [filename.split('/')[-3] for filename in protein_reports]  # samples names are the folder containing the reports folder 
        spectra_count = pd.DataFrame(columns = ['Main Accession'])
        for i in range(len(protein_reports)):
            report = pd.read_csv(protein_reports[i], sep = '\t', index_col = 0)
            if not uniprot_ids:
                blast = parse_blast(blast)
                blast['sseqid'] = [ide.split('|')[-1] for ide in blast.sseqid]
                report = pd.merge(report, blast[['qseqid','sseqid']], 
                                 left_on = 'Main Accession', right_on = 'qseqid')
                report = report[['sseqid', '#PSMs']]
            else:
                report = report[['Main Accession', '#PSMs']]
            report.columns = ['Main Accession', samples_names[i]]
            report = report.groupby('Main Accession')[samples_names[i]].sum().reset_index()
            spectra_count = pd.merge(spectra_count, report, on = 'Main Accession',
                                     how = 'outer')
        spectra_count[samples_names] = spectra_count[samples_names].fillna(value = 0).astype(int)
        spectra_count.to_csv(output, sep = '\t', index = False)
        
    '''
    MetaProteomics with MaxQuant
    '''
    
    '''
    input:
        output: name of standard parameters file to create
    output:
        a standard parameters file for MaxQuant named "output" will be created
    '''
    def create_mqpar(self, output):
        if os.path.isfile(output):                                              # the create file command will not create a new one if the file already exists
            os.remove(output)                                                   # even if that file already has non-default information in it, messing with the next commands
        run_command('maxquant ' + output + ' --create')
        
    '''
    input:
        mqpar: name of the mqpar.xml file to have its values changed
        fasta_database: name of the FASTA database for metaproteomics
        spectra_folder: name of the folder with RAW spectra files
        experiment_names: list of experiments, with one element for each file
    output:
        the "file" file will be updated with the new parameters
    '''
    def edit_maxquant_mqpar(self, mqpar, fasta_database, spectra_folder, 
                            experiment_names, threads = 1, file_type = 'raw'):
        print('Updating parameters file information.')
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(mqpar, parser)
        root = tree.getroot()
        root.find("fastaFiles/FastaFileInfo/fastaFilePath").text = fasta_database
        print('Fasta database = ' + fasta_database)
        root.find("separateLfq").text = 'True'
        root.find("numThreads").text = str(threads)
        print('Number of threads = ' + str(threads))
        for child in ['filePaths/string', 'experiments/string',
                                    'fractions/short', 'ptms/boolean', 
                                    'paramGroupIndices/int']:
            tree.xpath(child)[0].getparent().remove(tree.xpath(child)[0])       # the default params file of MaxQuant brings default arguments that can't be present
        filePaths = root.find("filePaths")
        experiments = root.find("experiments")
        fractions = root.find("fractions")
        ptms = root.find("ptms")
        paramGroupIndices = root.find("paramGroupIndices")
        files = sort_alphanumeric(glob.glob('{}/*.{}'.format(spectra_folder, file_type)))
        for i in range(len(files)):
            print('Adding file: ' + files[i])
            etree.SubElement(filePaths, 'string').text = files[i]
            print('Experiment = ' + experiment_names[i])
            etree.SubElement(experiments, 'string').text = experiment_names[i]
            etree.SubElement(fractions, 'short').text = '32767'
            etree.SubElement(ptms, 'boolean').text = 'True'
            etree.SubElement(paramGroupIndices, 'int').text = '0'
        root.find("parameterGroups/parameterGroup/lfqMode").text = '1'
        tree.write(mqpar, pretty_print=True)
        print('Parameters file is available at ' + mqpar)
        
    '''
    input:
        mqpar: name of MaxQuant parameters file
        spectra_folder: folder containing the spectra RAW files
        output_folder: folder where the "combined" folder will be placed
    output:
        
    '''
    def run_maxquant(self, mqpar, spectra_folder, output_folder):
        for directory in [spectra_folder + '/combined', output_folder + '/combined']:
            if os.path.isdir(directory):
                shutil.rmtree(directory, ignore_errors=True)
        run_command('maxquant ' + mqpar)        # TODO - get the shell messages from MaxQuant to appear
        os.rename(spectra_folder + '/combined', output_folder + '/maxquant_results')

    def compomics_workflow(self):
        self.verify_crap_db(self.crap_database)

        self.create_decoy_database(self.database, self.searchgui_exe)
        try:  # try/except - https://github.com/compomics/searchgui/issues/217
            self.generate_parameters_file(self.output + '/params.par',
                                          self.database.replace('.fasta', '_concatenated_target_decoy.fasta'),
                                          self.searchgui_exe)
        except:
            print('An illegal reflective access operation has occurred. But MOSCA can handle it.')
        self.peptide_spectrum_matching(self.spectra_folder, self.output,
                                       self.output + '/params.par',
                                       self.searchgui_exe, threads=self.threads)
        self.browse_identification_results(self.spectra_folder, self.output + '/params.par',
                                           self.output + '/searchgui_out.zip', self.output + '/ps_output.cpsx',
                                           self.peptide_shaker_exe)
        try:  # try/except - if no identifications are present, will throw an error
            self.generate_reports(self.output + '/ps_output.cpsx', self.output + '/reports',
                                  self.peptide_shaker_exe)
        except:
            print('No identifications?')

        self.spectra_counting((self.output + '/reports/' + self.experiment_name +
                               '_' + self.sample_name + '_' + self.replicate_number +
                               '_Default_Protein_Report.txt'), self.blast,
                              self.output + '/Spectra_counting.tsv')

    def maxquant_workflow(self):
        self.create_mqpar(self.output + '/mqpar.xml')
        self.edit_maxquant_mqpar(self.output + '/mqpar.xml', self.output + '/database.fasta',
                                 self.spectra_folder, self.experiment_names,
                                 threads=self.threads)
        self.run_maxquant(self.output + '/mqpar.xml', self.spectra_folder, self.output)
        
    def run(self):

        self.database_generation(self.output + '/database_uniprot_sequences.fasta', 
                                 self.crap_database, self.protease, 
                                 how = 'raw', blast = self.blast, faa = self.faa)
        
        if self.workflow == 'maxquant':
            self.maxquant_Workflow()
            
        elif self.workflow == 'compomics':
            self.compomics_workflow()
            
    def background_inputation(self, df):
        return df.fillna(value = df.min().min())
            
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
                df.loc[line][replicate] = df.loc[line][replicate].fillna(value = min_val)
    return df
                    
def background_inputation(df):
    return df.fillna(value = df.min())

if __name__ == '__main__':
    MetaproteomicsAnalyser().run()