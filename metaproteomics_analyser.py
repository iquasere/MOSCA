# -*- coding: utf-8 -*-
"""
MOSCA's MetaProteomics class for performing MetaProteomics Analysis

By João Sequeira

Jul 2018
"""

from progressbar import ProgressBar
from mosca_tools import MoscaTools
from diamond import DIAMOND
import pandas as pd
from lxml import etree
import os, pathlib, shlex, subprocess, glob

mtools = MoscaTools()

class MetaProteomicsAnalyser:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    '''   
    input: 
        faa: fasta file from FragGeneScan
        output: output folder for protease file and final database
        crap: filename of cRAP database
        protease: protease name used
    output: 
        a FASTA file named "output + /database.fasta" will be created,
        containing the sequences from metagenomics, from cRAP database and of the
        protease
    '''
    def database_generation(self, faa, output, crap, protease = 'trypsin'):      # SearchCLI doesn't allow for repeated IDs, so for now will not use them # TODO - find a workaround
        handler = open(output + '/db.fasta', 'w')
        faa = mtools.parse_fasta(faa)
        for k,v in faa.items():
            if '*' not in v:
                handler.write('>' + k + '\n' + v + '\n')
        '''
        blast = DIAMOND(out = blast).parse_result()
        faa = pd.DataFrame.from_dict(faa, orient = 'index')
        faa.columns = ['sequence']
        database = pd.merge(blast, faa, left_on = 'qseqid', right_index = True)
        database = database[['sseqid', 'sequence']]
        database['sseqid'] = [ide.split('|')[-1] for ide in database.sseqid]
        pbar = ProgressBar()
        print('Generating new database to ' + output + '/database.fasta')
        for ide in pbar(database.index):
            handler.write('>' + database.loc[ide]['sseqid'] + '\n' + database.loc[ide]['sequence'] + '\n')
        '''
        handler.close()
        if protease == 'trypsin':                                                                   #pig trypsin is not the only protease used # TODO - will have to include more
            if not os.path.isfile(output + '/P00761.fasta'):
                print('Trypsin file not found. Will be downloaded from UniProt.')
                mtools.run_command('wget https://www.uniprot.org/uniprot/P00761.fasta -P ' + output)
            protease = output + '/P00761.fasta'
        mtools.run_command('cat ' + ' '.join([output + '/db.fasta', crap, protease]), 
                           file = output + '/database.fasta', mode = 'w')
        os.remove(output + '/db.fasta')
    
    '''   
    input: 
        crap_folder: folder where crap.fasta should be
    output: 
        confirmation of cRAP database presence in the folder, or download of the 
        database if absent
    '''    
    def verify_crap_db(self, crap_folder):
        if os.path.isfile(crap_folder + '/crap.fasta'):
            print('cRAP database exists at ' + crap_folder + '/crap.fasta.')
        else:
            print('cRAP database not found at ' + crap_folder + '/crap.fasta. Downloading cRAP database.')
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
        bashCommand = ('java -cp ' + os.path.expanduser('~/SearchGUI-3.3.3/SearchGUI-3.3.3.jar') + 
                       ' eu.isas.searchgui.cmd.FastaCLI -in ' + database + ' -decoy')
        mtools.run_command(bashCommand)
        
    '''   
    input: 
        output: name of parameters file
        database: name of FASTA decoy database
    output:
        a parameters file will be produced for SearchCLI and/or PeptideShakerCLI
    ''' 
    def generate_parameters_file(self, output, database):
        bashCommand = ('java -cp ' + os.path.expanduser('~/SearchGUI-3.3.3/SearchGUI-3.3.3.jar') + 
                       ' eu.isas.searchgui.cmd.IdentificationParametersCLI -out ' + 
                       output + ' -db ' + database + ' -prec_tol 10 -frag_tol 0.02' +
                       ' -enzyme Trypsin -fixed_mods "Carbamidomethylation of C"' +
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
        bashCommand = ('java -cp ' + os.path.expanduser('~/SearchGUI-3.3.3/SearchGUI-3.3.3.jar') + 
                        ' eu.isas.searchgui.cmd.SearchCLI -spectrum_files ' + spectra_folder + 
                        ' -output_folder ' + output + ' -id_params ' + parameters_file +
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
        bashCommand = ('java -cp ' + os.path.expanduser('~/PeptideShaker-1.16.31/PeptideShaker-1.16.31.jar') + 
                       ' eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files ' +
                       spectra_folder + ' -experiment ' + experiment_name + ' -sample ' + 
                       sample_name + ' -replicate ' + replicate_number + ' -identification_files ' + 
                       searchcli_output + ' -out ' + peptideshaker_output)
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
        bashCommand = ('java -cp ' + os.path.expanduser('~/PeptideShaker-1.16.31/PeptideShaker-1.16.31.jar') + 
                       ' eu.isas.peptideshaker.cmd.ReportCLI -in ' + peptideshaker_output + 
                       ' -out_reports ' + reports_folder + ' -reports ' + ','.join(reports_list))
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
        
    def run(self):
        self.verify_crap_db(self.crap_folder)
        self.database_generation(self.faa, self.output, self.crap_folder + '/crap.fasta', self.protease)
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
        mtools.run_command('mono ' + os.path.expanduser('~/MaxQuant/bin/MaxQuantCmd.exe ' +
                                                        output + ' --create'))
        
    '''
    input:
        file: name of the mqpar.xml file to have its values changed
        fasta_database: name of the FASTA database for metaproteomics
        spectra_folder: name of the folder with RAW spectra files
        mg_name: name of the experiment
    output:
        the "file" file will be updated with the new parameters
    '''
    def edit_maxquant_mqpar(self, mqpar, fasta_database, spectra_folder, mg_name, 
                            threads = 1):
        print('Updating parameters file information.')
        parser = etree.XMLParser(remove_blank_text=True)
        tree = etree.parse(mqpar, parser)
        root = tree.getroot()
        root.find("fastaFiles/FastaFileInfo/fastaFilePath").text = fasta_database
        print('Fasta database = ' + fasta_database)
        root.find("numThreads").text = str(threads)
        print('Number of threads = ' + str(threads))
        files = mtools.sort_alphanumeric(glob.glob(spectra_folder + '/*.RAW'))
        print('Spectra folder = ' + spectra_folder)
        for child in ['filePaths/string', 'experiments/string',
                                'fractions/short', 'ptms/boolean', 
                                'paramGroupIndices/int']:
            tree.xpath(child)[0].getparent().remove(tree.xpath(child)[0])       # the default params file of MaxQuant brings default arguments that can't be present
        filePaths = root.find("filePaths")
        experiments = root.find("experiments")
        fractions = root.find("fractions")
        ptms = root.find("ptms")
        paramGroupIndices = root.find("paramGroupIndices")
        for file in files:
            print('Adding information for spectra file ' + file)
            etree.SubElement(filePaths, 'string').text = file
            etree.SubElement(experiments, 'string').text = mg_name
            etree.SubElement(fractions, 'short').text = '32767'
            etree.SubElement(ptms, 'boolean').text = 'False'
            etree.SubElement(paramGroupIndices, 'int').text = '0'
        tree.write(mqpar, pretty_print=True)
        print('Parameters file is available at ' + mqpar)
        
    '''
    input:
        mqpar: name of MaxQuant parameters file
    output:
        
    '''
    def run_maxquant(self, mqpar):
        mtools.run_command('mono ' + os.path.expanduser('~/MaxQuant/bin/MaxQuantCmd.exe ' +         # TODO - get the shell messages from MaxQuant to appear
                                                        mqpar))
        
if __name__ == '__main__':
    correspondence = {'PAL6_S2_L001':'DNA7'}  # 'EST6_S1_L001':'DNA4','OL6_S3_L001':'DNA5','OLDES6_S4_L001':'DNA6',
    
    mper = MetaProteomicsAnalyser()
    
    for experiment in correspondence.keys():
        path = pathlib.Path('MGMP/MetaProteomics/' + experiment).mkdir(parents=True, exist_ok=True)
        mper.create_mqpar('/HDDStorage/jsequeira/Thesis/MGMP/MetaProteomics/' + experiment + '/mqpar.xml')
        mper.edit_maxquant_mqpar('/HDDStorage/jsequeira/Thesis/MGMP/MetaProteomics/' + experiment + '/mqpar.xml',
                                 '/HDDStorage/jsequeira/Thesis/MGMP/MetaProteomics/' + experiment + '/database.fasta',
                                 '/HDDStorage/jsequeira/Thesis/Datasets/Proteic/' + correspondence[experiment],
                                 experiment, threads = 6)
        mper.run_maxquant('/HDDStorage/jsequeira/Thesis/MGMP/MetaProteomics/' + experiment + '/mqpar.xml')
    
    '''
    for sample in correspondence.keys():
        
        pathlib.Path('MGMP/' + sample).mkdir(parents=True, exist_ok=True)
        
        analyser = MetaProteomicsAnalyser(faa = 'MGMP/Annotation/' + sample + '/fgs.faa',
                                          blast = 'MGMP/Annotation/' + sample + '/aligned.blast',
                                          crap_folder = '/HDDStorage/jsequeira/Thesis/MetaProteomics',
                                          output = 'MGMP/Analysis/' + sample,
                                          protease = 'trypsin',
                                          spectra_folder = 'Datasets/Proteic/' + correspondence[sample],
                                          experiment_name = 'MGMP',
                                          sample_name = sample,
                                          replicate_number = '1')
        
        analyser.run()
    '''
    '''
    import glob
    
    files = glob.glob('MGMP/Analysis/*/Spectra_counting.tsv')
    
    read_spectra = pd.DataFrame()
    
    readspectra = pd.read_csv(files[0], sep = '\t', header = None, index_col = 0)
    print(readspectra.head())
    readspectra.columns = [files[0].split('/')[-2]]
    
    for file in files[1:]:
        df = pd.read_csv(file, sep = '\t', header = None, index_col = 0)
        df.columns = [file.split('/')[-2]]
        readspectra = pd.merge(readspectra, df, left_index = True, right_index = True, how = 'outer')
    readspectra = readspectra.fillna(value = 0)
    readspectra.index.name = 'geneid'
    readspectra.to_csv('MGMP/Analysis/Spectra_counting.tsv',sep='\t')
    '''