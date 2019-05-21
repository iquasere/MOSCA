# -*- coding: utf-8 -*-
"""
MOSCA's MetaProteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

from progressbar import ProgressBar
from mosca_tools import MoscaTools
from diamond import DIAMOND
from annotation import Annotater
import pandas as pd
from lxml import etree
import os, pathlib, shlex, subprocess, glob, shutil

mtools = MoscaTools()

class MetaProteomicsAnalyser:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
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
    def database_generation(self, output, crap_database, protease = 'trypsin',
                            how = 'raw', faa = None, blast = None):
        temp = '/'.join(output.split('/')[:-1]) + '/temporary.fasta'
        print('Generating new database to ' + output)
        if how == 'uniprot_sequences':
            annotater = Annotater()
            annotater.recursive_uniprot_fasta(output, blast = blast)
            
        elif how == 'raw':
            database = list()
            faa = mtools.parse_fasta(faa)
            print('Removing sequences with unknown bases (*).')
            pbar = ProgressBar()
            for k,v in pbar(faa.items()):
                if '*' not in v:
                    database.append([k,v])
                    
        elif how == 'uniprot_ids':
            faa = mtools.parse_fasta(faa)
            faa = pd.DataFrame.from_dict(faa, orient='index')
            faa.columns = ['sequence']
            faa = pd.merge(blast, faa, left_on = 'qseqid', right_index = True)
            with open(temp, 'w') as file:
                print('Writing database with UniProt IDs to ' + temp)
                pbar = ProgressBar()
                for i in pbar(range(len(faa))):
                    file.write('>' + faa.iloc[i]['sseqid'] + '\n' + faa.iloc[i]['sequence'] + '\n')
        if protease == 'trypsin':                                                                   # TODO - pig trypsin is not the only protease used, will have to include more - see from searchgui
            if not os.path.isfile(output + '/P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                mtools.run_command('wget https://www.uniprot.org/uniprot/P00761.fasta' + 
                                   ' -P ' + '/'.join(output.split('/')[:-1]))
            protease = output + '/P00761.fasta'
        mtools.run_command('cat ' + ' '.join([crap_database, protease]), 
                           file = output, mode = 'a')

    '''   
    input: 
        crap_folder: folder where crap.fasta should be
    output: 
        confirmation of cRAP database presence in the folder, or download of the 
        database if absent
    '''    
    def verify_crap_db(self, crap_database):
        if os.path.isfile(crap_database):
            print('cRAP database exists at ' + crap_database)
        else:
            print('cRAP database not found at ' + crap_database + '. Downloading cRAP database.')
            crap_folder = '/'.join(crap_database.split('/')[:-1])
            mtools.run_command('wget ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta -P ' + crap_folder)
            
    '''   
    input:
        protein_fasta: fasta file with proteins from MG analysis, plus trypsin 
        and cRAP sequences
    output:
        a FASTA file named "database + _concatenated_target_decoy.fasta" will be 
        created, containing interleaved original and decoy sequences
    ''' 
    def create_decoy_database(self, database):
        bashCommand = ('searchgui eu.isas.searchgui.cmd.FastaCLI -in ' 
                       + database + ' -decoy')
        mtools.run_command(bashCommand)
        
    '''   
    input: 
        output: name of parameters file
        database: name of FASTA decoy database
    output:
        a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
    ''' 
    def generate_parameters_file(self, output, database):
        bashCommand = ('searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI ' + 
                       '-out ' + output + ' -db ' + database + ' -prec_tol 10 ' + 
                       '-frag_tol 0.02 -enzyme Trypsin -fixed_mods "Carbamidomethylation of C"' +
                       ' -variable_mods "Oxidation of M" -mc 2')
        print(bashCommand)
        bashCommand_correct = shlex.split(bashCommand)                              #the usual way with MoscaTools is not working here, dunno why
        process = subprocess.Popen(bashCommand_correct, stdout=subprocess.PIPE)
        output, error = process.communicate()
        #mtools.run_command(bashCommand)
        
    '''   
    input: 
        spectra_folder: folder containing the raw spectra files
        output: folder to output results
        parameters_file: parameters filename 
    output:
        a "searchgui_out.zip" file will be created in the output folder
    ''' 
    def peptide_spectrum_matching(self, spectra_folder, output, parameters_file,
                                  search_engines = ['xtandem', 'myrimatch', 'msgf']):
        bashCommand = ('searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files '
                       + spectra_folder + ' -output_folder ' + output + ' -id_params ' 
                       + parameters_file +
                       ''.join([' -' + engine + ' 1' for engine in search_engines]))
        mtools.run_command(bashCommand)
        
    '''   
    input: 
        spectra_folder: folder containing the raw spectra files
        experiment_name: name of experiment
        sample_name: name of sample
        replicate_number: number of replicate (STRING)
        output: folder to output results
        parameters_file: parameters filename
        searchcli_output: searchcli output filename
        peptideshaker_output: peptideshaker output filename
    output:
        a file will be outputed with the validation of the PSMs and protein
        identifications
    ''' 
    def browse_identification_results(self, spectra_folder, experiment_name, 
                                      sample_name, replicate_number, parameters_file,
                                      searchcli_output, peptideshaker_output):
        bashCommand = ('peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI ' +
                       '-spectrum_files ' + spectra_folder + ' -experiment ' + 
                       experiment_name + ' -sample ' + sample_name + ' -replicate ' + 
                       replicate_number + ' -identification_files ' + searchcli_output + 
                       ' -out ' + peptideshaker_output)
        mtools.run_command(bashCommand)
       
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
    def generate_reports(self, peptideshaker_output, reports_folder, 
                         reports_list = [str(n) for n in range(12)]):
        print('Created ' + reports_folder)
        pathlib.Path(reports_folder).mkdir(parents=True, exist_ok=True)                 # creates folder for reports
        bashCommand = ('peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in ' + 
                       peptideshaker_output + ' -out_reports ' + reports_folder + 
                       ' -reports ' + ','.join(reports_list))
        mtools.run_command(bashCommand)
        
    '''   
    input: 
        protein_report: name of file containing protein report from PeptideShaker
        blast: name of blast file with UniProt annotations
        output: name of file to output
    output:
        A tab separated spectra count file named "output", with UniProt IDs
        and corresponding spectra quantification
    '''
    def spectra_counting(self, protein_report, blast, output):
        report = pd.read_csv(protein_report, sep = '\t', index_col = 0)
        blast = DIAMOND(out = blast).parse_result()
        blast['sseqid'] = [ide.split('|')[-1] for ide in blast.sseqid]
        report = pd.merge(report, blast[['qseqid','sseqid']], 
                         left_on = 'Main Accession', right_on = 'qseqid')
        report = report[['sseqid', '#Peptides']]
        report = report.groupby('sseqid')['#Peptides'].sum().reset_index()
        report.to_csv(output, sep = '\t', index = False, header = None)
        
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
        mtools.run_command('maxquant' + output + ' --create')
        
    '''
    input:
        file: name of the mqpar.xml file to have its values changed
        fasta_database: name of the FASTA database for metaproteomics
        spectra_folder: name of the folder with RAW spectra files
        experiment_name: name of the experiment
    output:
        the "file" file will be updated with the new parameters
    '''
    def edit_maxquant_mqpar(self, mqpar, fasta_database, spectra_folder, 
                            experiment_name, threads = 1):
        print('Updating parameters file information.')
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(mqpar, parser)
        root = tree.getroot()
        root.find("fastaFiles/FastaFileInfo/fastaFilePath").text = fasta_database
        print('Fasta database = ' + fasta_database)
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
        print('Experiment = ' + experiment_name)
        files = mtools.sort_alphanumeric(glob.glob(spectra_folder + '/*.RAW'))
        for file in files:
            print('Adding file: ' + file)
            etree.SubElement(filePaths, 'string').text = file
            etree.SubElement(experiments, 'string').text = experiment_name
            etree.SubElement(fractions, 'short').text = '32767'
            etree.SubElement(ptms, 'boolean').text = 'False'
            etree.SubElement(paramGroupIndices, 'int').text = '0'
        '''
        for experiment in experiments_files.index:
            print('Experiment = ' + experiment)
            files = mtools.sort_alphanumeric(experiments_files.loc[experiment]['files'])
            for file in files:
                print(file)
                etree.SubElement(filePaths, 'string').text = file
                etree.SubElement(experiments, 'string').text = experiment
                etree.SubElement(fractions, 'short').text = '32767'
                etree.SubElement(ptms, 'boolean').text = 'False'
                etree.SubElement(paramGroupIndices, 'int').text = '0'
        '''
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
        mtools.run_command('maxquant ' + mqpar)        # TODO - get the shell messages from MaxQuant to appear
        #os.rename(spectra_folder + '/combined', output_folder + '/maxquant_results')
        
    def run(self):
        
        self.verify_crap_db(self.crap_database)
        self.database_generation(self.output + '/database_uniprot_sequences.fasta', 
                                 self.crap_database, self.protease, 
                                 how = 'raw', blast = self.blast, faa = self.faa)
        
        if self.workflow == 'maxquant':
            self.create_mqpar(self.output + '/mqpar.xml')
            self.edit_maxquant_mqpar(self.output + '/mqpar.xml', self.output + '/database.fasta',
                                     self.spectra_folder, self.experiment_name,
                                     threads = 6)
            self.run_maxquant(self.output + '/mqpar.xml', self.spectra_folder, self.output)
            
        elif self.workflow == 'compomics':
            self.create_decoy_database(self.output + '/database.fasta')
            self.generate_parameters_file(self.output + '/params.par', 
                                          self.output + '/database_concatenated_target_decoy.fasta')
            self.peptide_spectrum_matching(self.spectra_folder, self.output, 
                                           self.output + '/params.par')
            self.browse_identification_results(self.spectra_folder, self.experiment_name, 
                        self.sample_name, self.replicate_number, self.output + '/params.par', 
                        self.output + '/searchgui_out.zip', self.output + '/ps_output.cpsx')
            self.generate_reports(self.output + '/ps_output.cpsx', self.output + '/reports')
            self.spectra_counting((self.output + '/reports/' + self.experiment_name + 
                                  '_' + self.sample_name + '_' +  self.replicate_number + 
                                  '_Default_Protein_Report.txt'), self.blast,
                                  self.output + '/Spectra_counting.tsv')