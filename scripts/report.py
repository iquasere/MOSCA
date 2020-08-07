# -*- coding: utf-8 -*-
'''
General report construction and export

By João Sequeira

Oct 2019
'''

from mosca_tools import MoscaTools
import pandas as pd, numpy as np, os, glob

mtools = MoscaTools()

class Reporter:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        
        self.fastq_columns = open('{}/fastqc_columns.txt'.format(self.lists_dir)).read().split('\n')
        self.metaquast_columns = open('{}/metaquast_columns.txt'.format(self.lists_dir)).read().split('\n')
        
    '''
    Input:
        output: str - filename to write tools and respective versions
    Output:
        a file named [output] will be written with information concerning
        the softwares used by MOSCA and respective versions
    '''
    def write_technical_report(self, output):                                   # TODO - add proteomics software that cannot be installed with conda
        conda_list = mtools.run_pipe_command('conda list', output = 'PIPE').split('\n')[2:]
        lines = [line.split() for line in conda_list]
        lines[0] = lines[0][1:]
        df = pd.DataFrame(lines, columns = lines.pop(0)).set_index('Name')
        tools = open('{}/tools_for_versions.txt'.format(self.lists_dir)).read().split('\n')
        present_tools = [tool for tool in tools if tool in df.index]
        open(output, 'w').write(df.loc[present_tools][['Version']].to_string(
            justify = 'left'))
        
    def initialize_report(self):
        print('Initializing Report')
        self.report = pd.DataFrame(columns = open(
            '{}/reporter_columns.txt'.format(self.lists_dir)).read().split('\n'))
    
    '''
    Input: 
        name: str - name of MG/MT sample
        prefix: str - [Initial quality assessment], [Before quality trimming] or
        [After quality trimming] - the step that concerns the FastQC report
        performed_rrna_removal: bool - True if rRNA removal was performed, False
        otherwise. Should be False for MG, and True for MT
    Output:
        self.report will be updated with information from the report, on the line
        named 'name', and the columns that start with 'prefix'
    '''
    def info_from_fastqc(self, output_dir, name, prefix, prefix2terms):
        reports = [mtools.parse_fastqc_report(
                '{}/Preprocess/FastQC/{}{}_{}_fastqc/fastqc_data.txt'.format(output_dir, 
                 prefix2terms[prefix][0], name, prefix2terms[prefix][i])) for i in [1, 2]]
        for column in self.fastq_columns:
            for i in range(2):
                if column not in reports[i].keys():
                    reports[i][column] = ('Not available', None)                # only the not available matters. And nothing else matters!...
            if reports[0][column][0] == reports[1][column][0]:
                self.report.loc[name, '{} {}'.format(prefix, column)] = reports[0][column][0]
            else:
                self.report.loc[name, '{} {}'.format(prefix, column)] = (
                        '{} (forward) {} (reverse)'.format(
                                reports[0][column][0], reports[1][column][0]))
                        
    def info_from_preprocessing(self, output_dir, name, input_file, performed_rrna_removal = False):
        print('Retrieving preprocessing information for sample: ' + name)
        if name not in self.report.index:
            self.report = self.report.append(pd.Series(name = name))
        self.report.loc[name] = self.report.loc[name].fillna(value = '')
        
        adapter_files = open('{}/Preprocess/Trimmomatic/{}_adapters.txt'.format(output_dir, name)).read().split('\n')
        if len(adapter_files[0]) > 0:  
            adapter = adapter_files[0].split('/')[-1].split('.fa')[0]
        else:
           adapter_files = list(); adapter = None
           
        # For each preprocessing step, a tuple of (prefix, suffix for forward, suffix for reverse)
        prefix2terms = {'[Initial quality assessment]': ('', 'R1', 'R2'),
             '[Before quality trimming]': (('', 'forward', 'reverse') 
             if performed_rrna_removal else ('', adapter + '_forward_paired', 
                                               adapter + '_reverse_paired') 
             if adapter is not None else ('', 'R1', 'R2')), '[After quality trimming]': (
                     'quality_trimmed_', 'forward_paired', 'reverse_paired')}
        '''
        # Initial assessment
        self.report.loc[name, '[Initial quality assessment] # of initial reads'] = mtools.count_on_file(
            '@', input_file, compressed = True if input_file.endswith('.gz') else False)
        self.info_from_fastqc(output_dir, name, '[Initial quality assessment]', prefix2terms)
        
        # After adapter removal
        try:
            if len(adapter_files) > 0:
                self.report.loc[name, '[Adapter removal] adapter files'] = ', '.join(set(adapter_files))
                self.report.loc[name, '[Adapter removal] # of reads remaining'] = (
                    mtools.count_on_file('@', '{}/Preprocess/Trimmomatic/{}_{}_forward_paired.fq'.format(
                        output_dir, name, adapter)))
            else:
                self.report.loc[name, '[Adapter removal] adapter files'] = 'None'
                self.report.loc[name, '[Adapter removal] # of reads remaining'] = (
                    self.report.loc[name]['[Initial quality assessment] # of initial reads'])
        except:
            print('Failed at adapter removal!')
            self.report.to_csv('{}/report.tsv'.format(output_dir),sep='\t')
            
        # rRNA removal
        if performed_rrna_removal:
            self.report.loc[name, '[rRNA removal] # of reads remaining'] = (
            mtools.count_on_file('@', '{}/Preprocess/SortMeRNA/{}_forward.fastq'.format(
                    output_dir, name)))
        else:
            self.report.loc[name, '[rRNA removal] # of reads remaining'] = (
                self.report.loc[name]['[Adapter removal] # of reads remaining'])
        
        # Quality trimming
        self.info_from_fastqc(output_dir, name, '[Before quality trimming]', prefix2terms)
        '''
        self.report.loc[name, '[Quality trimming] Parameters'] = '; '.join([file for file in 
                set(open('{}/Preprocess/Trimmomatic/{}_quality_params.txt'.format(
                output_dir, name)).read().split('\n')) if len(file) > 2])                                # TODO - because '' must be interfering, try to cut the problem at the root before troubles
        self.report.loc[name, '[Quality trimming] # of reads remaining'] = (
            mtools.count_on_file('@', '{}/Preprocess/Trimmomatic/quality_trimmed_{}_forward_paired.fq'.format(
                    output_dir, name)))
        self.info_from_fastqc(output_dir, name, '[After quality trimming]', prefix2terms)
        self.report.to_csv('{}/report1.tsv'.format(output_dir), sep='\t')
        
    def set_samples(self, mg_name2sample, mg2mt):
        mg_name2sample = [[mg_name, sample] for mg_name, sample in mg_name2sample.items()]
        mtname2sample = list()
        for pair in mg_name2sample:
            for mt_name in mg2mt[pair[0]]:
                mtname2sample.append([mt_name, pair[1]])
        name2sample = mg_name2sample + mtname2sample
        name2sample = pd.DataFrame(name2sample, columns = ['Name', 'Sample'])
        print(self.report)
        self.report = pd.merge(name2sample, self.report, left_on = 'Name', 
                               right_index = True, how = 'outer')
        print(self.report)
    
    def info_from_assembly(self, output_dir, sample):
        print('Retrieving assembly information for sample ' + sample)
        qc_report = pd.read_csv('{}/Assembly/{}/quality_control/report.tsv'.format(
                output_dir, sample), sep = '\t', index_col = 0).transpose()
        qc_report.index = [sample]
        for col in qc_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Assembly] ' + col] = (
                    qc_report.loc[sample][col])
        self.report.to_csv('{}/report2.tsv'.format(output_dir),sep='\t')
        
    def info_from_annotation(self, output_dir, sample):
        print('Retrieving annotation information for sample ' + sample)
        sample_report = dict()
        sample_report['# of proteins detected'] = (
            mtools.count_on_file('>', '{}/Annotation/{}/fgs.faa'.format(output_dir, sample)))
        
        sample_report['# of proteins annotated (DIAMOND)'] = (
            len(set(mtools.parse_blast('{}/Annotation/{}/aligned.blast'.format(
                output_dir, sample))['qseqid'])) - mtools.count_on_file('*', 
                    '{}/Annotation/{}/aligned.blast'.format(output_dir, sample)))
        
        sample_report['# of proteins annotated (reCOGnizer)'] = (
            len(set(mtools.parse_blast('{}/Annotation/{}/cdd_aligned.blast'.format(
                    output_dir, sample))['qseqid'])))
        
        sample_report = pd.DataFrame.from_dict(sample_report, orient = 'index').transpose()
        sample_report.index = [sample]
        
        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Annotation] ' + col] = (
                    sample_report.loc[sample][col])
        self.report.to_csv('{}/report3.tsv'.format(output_dir),sep='\t')
    
    def info_from_binning(self, output_dir, sample):
        sample_report = dict()
        sample_report['# of bins'] = len(glob.glob(
            '{0}/Binning/{1}/{1}.*.fasta'.format(output_dir, sample)))
        checkm = pd.read_csv('{}/Binning/{}/CheckM_results/output.tab'.format(
            output_dir, sample), sep = '\t')
        sample_report['# of high-quality drafts'] = ((checkm['Completeness'] > 90) 
            & (checkm['Contamination'] < 5)).sum()
        sample_report['# of medium-quality drafts'] = ((checkm['Completeness'] < 90)
            & (checkm['Completeness'] > 50) & (checkm['Contamination'] < 10)).sum()
        sample_report['# of low-quality drafts'] = ((checkm['Completeness'] < 50) 
            & (checkm['Contamination'] < 10)).sum()
        sample_report = pd.DataFrame.from_dict(sample_report, orient = 'index').transpose()
        sample_report.index = [sample]
        
        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Binning] ' + col] = (
                    sample_report.loc[sample][col])
        self.report.to_csv('{}/report4.tsv'.format(output_dir),sep='\t')
        
    def info_from_alignment(self, output_dir, mt_name):
        self.report.set_index('Name', inplace = True)
        self.report.loc[mt_name, '[Metatranscriptomics] # of reads aligned'] = (
            mtools.run_pipe_command("""samtools view {}/Metatranscriptomics/{}.sam | 
                cut -f 3 | sort | uniq -c | awk \'{{printf(\"%s\\t%s\\n\", $2, $1)}}\' | 
                awk '{{sum+=$2}} END {{print sum}}\'""".format(output_dir, mt_name), 
                output = 'PIPE').rstrip('\n'))
        self.report.reset_index(inplace = True)
        self.report.to_csv('{}/report5.tsv'.format(output_dir),sep='\t')
            
    def info_from_differential_expression(self, output_dir, sample, cutoff = 0.01):
        de_results = pd.read_csv('{}/Metatranscriptomics/condition_treated_results.csv'.format(
            output_dir), index_col = 0)
        self.report.loc[self.report['Sample'] == sample, '[Gene expression] # of differentially expressed proteins'] = (
                    (de_results['padj'] < cutoff).sum())
        self.report.to_csv('{}/report6.tsv'.format(output_dir),sep='\t')
