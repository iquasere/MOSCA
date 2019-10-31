# -*- coding: utf-8 -*-
'''
General report construction and export

By João Sequeira

Oct 2019
'''

from mosca_tools import MoscaTools
import pandas as pd, numpy as np, os

mtools = MoscaTools()

class Reporter:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        
        self.fastq_columns = ['Per base sequence quality', 'Per tile sequence quality',
                 'Per sequence quality scores', 'Per base sequence content',
                 'Per sequence GC content', 'Per base N content',
                 'Sequence Length Distribution', 'Sequence Duplication Levels', 
                 'Overrepresented sequences', 'Adapter Content']
        
        self.metaquast_columns = ['# contigs (>= 0 bp)', '# contigs (>= 1000 bp)', 
                     '# contigs (>= 5000 bp)', '# contigs (>= 10000 bp)',
                     '# contigs (>= 25000 bp)', '# contigs (>= 50000 bp)',
                     'Total length (>= 0 bp)', 'Total length (>= 1000 bp)',
                     'Total length (>= 5000 bp)', 'Total length (>= 10000 bp)',
                     'Total length (>= 25000 bp)', 'Total length (>= 50000 bp)',
                     '# contigs', 'Largest contig', 'Total length', 
                     'Reference length', 'N50', 'N75', 'L50', 'L75',
                     '# misassemblies', '# misassembled contigs', 
                     'Misassembled contigs length', '# local misassemblies',
                     '# unaligned mis. contigs', '# unaligned contigs',
                     'Unaligned length', 'Genome fraction (%)',
                     'Duplication ratio', "# N's per 100 kbp",
                     '# mismatches per 100 kbp', '# indels per 100 kbp',
                     'Largest alignment', 'Total aligned length', 'Reads aligned (%)']
        
        self.taxonomic_columns = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                     'genus', 'species']
        
        self.report = self.initialize_report(self.metatranscriptomics)
        

    def initialize_report(self, metatranscriptomics):
        return pd.DataFrame(columns = (['Initial quality assessment] # of initial reads'] +
        ['[Initial quality assessment] ' + col for col in self.fastq_columns] + 
        ['[Adapter removal] adapter files', '[Adapter removal] # of reads remaining', 
        'Adapter removal] % of reads removed', '[rRNA removal] # of reads remaining',
        '[rRNA removal] % of reads removed'] + 
        ['[Before quality trimming] ' + col for col in self.fastq_columns] +
        ['[Quality trimming] # of reads remaining', '[Quality trimming] % of reads removed'] + 
        ['[After quality trimming] ' + col for col in self.fastq_columns] + 
        ['[Assembly] ' + col for col in self.metaquast_columns] +
        ['[Annotation] ' + col for col in (['# of proteins detected', 
         '# of proteins annotated (DIAMOND)', '% of proteins annotated (DIAMOND)',
         '# of proteins annotated (PSI-BLAST)', '% of proteins annotated (PSI-BLAST)'] +
         ['Main taxa identified (' + col + ')' for col in self.taxonomic_columns] + 
         ['Main functional categories identified'])] +
        ['[Binning] ' + col for col in ['# of bins', '# of useful bins',
         'average completeness', 'average contamination']] + 
        (['[Metatranscriptomics] ' + col for col in (['# of reads aligned', 
          '% of reads aligned', '# of differentially expressed proteins'] + 
         ['Main taxa identified (' + col + ')' for col in self.taxonomic_columns] +
         ['Main functional categories identified'])] if metatranscriptomics else 
         ['[Metaproteomics] ' + col for col in (
                 ['# of differentially expressed proteins'] + 
         ['Main taxa identified (' + col + ')' for col in self.taxonomic_columns] +
         ['Main functional categories identified'])])))
    
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
    # TODO - for prefix in prefix2terms.keys(), do as stated under this
    def info_from_fastqc(self, output_dir, name, prefix, prefix2terms):
        reports = [mtools.parse_fastqc_report(
                '{}/Preprocess/FastQC/{}{}_{}_fastqc/fastqc_data.txt'.format(output_dir, 
                 prefix2terms[prefix][0], name, prefix2terms[prefix][i])) for i in [1, 2]]
        for column in self.fastq_columns:
            if reports[0][column][0] == reports[1][column][0]:
                self.report.at[name, '{} {}'.format(prefix, column)] = reports[0][column][0]
            else:
                self.report.at[name, '{} {}'.format(prefix, column)] = (
                        '{} (forward) {} (reverse)'.format(
                                reports[0][column][0], reports[1][column][0]))
                        
    def info_from_preprocessing(self, output_dir, name, input_file, performed_rrna_removal = False):
        self.report = self.report.append(pd.Series(name = name))
        self.report.loc[name] = self.report.loc[name].fillna(value = '')

        adapter_files = open('{}/Preprocess/Trimmomatic/{}_adapters.txt'.format(output_dir, name)).read().split('\n')
        if adapter_files[0] == ['']:  
            adapter = adapter_files[0].split('/')[-1]
        else:
           adapter_files = list(); adapter = None
        
        # For each preprocessing step, a tuple of (prefix, suffix for forward, suffix for reverse)
        prefix2terms = {'[Initial quality assessment]': ('', 'R1', 'R2'),
             '[Before quality trimming]': (('', 'forward', 'reverse') 
             if self.metatranscriptomics else ('after_adapter_removal_', '_' + 
                adapter + '_forward_paired', '_' + adapter + '_reverse_paired') 
             if adapter is not None else ('', 'R1', 'R2')), '[After quality trimming]': (
                     'quality_trimmed_', 'forward_paired', 'reverse_paired')}
        print('Initial assessment')
        # Initial assessment
        self.report.at[name, '[Initial quality assessment] # of initial reads'] = mtools.count_on_file(
            '@', input_file, compressed = True if input_file.endswith('.gz') else False)
        self.info_from_fastqc(output_dir, name, '[Initial quality assessment]', prefix2terms)
        print('After adapter removal')
        # After adapter removal
        if len(adapter_files) > 0:
            self.report.at[name, '[Adapter removal] adapter files'] = ', '.join(adapter_files)
            self.report.at[name, '[Adapter removal] # of reads remaining'] = (
                mtools.count_on_file('@', '{}/Preprocess/Trimmomatic/after_adapter_removal_{}_{}_forward_paired.fq'.format(
                    output_dir, name, adapter)))
            self.report.at[name, '[Adapter removal] % of reads removed'] = (
                    (self.report.loc[name]['# of initial reads'] - 
                     self.report.loc[name]['[Adapter removal] # of reads remaining']) /
                    self.report.loc[name]['# of initial reads'] * 100)
        else:
            self.report.at[name, '[Adapter removal] adapter files'] = 'None'
            self.report.at[name, '[Adapter removal] # of reads remaining'] = (
                self.report.loc[name]['[Initial quality assessment] # of initial reads'])
            self.report.at[name, '[Adapter removal] % of reads removed'] = 0.0
        print('rRNA removal')
        # rRNA removal
        if performed_rrna_removal:
            self.report.at[name, '[rRNA removal] # of reads remaining'] = (
            mtools.count_on_file('@', '{}/Preprocess/SortMeRNA/{}_forward.fastq'.format(
                    output_dir, name)))
            self.report.at[name, '[rRNA removal] % of reads removed'] = (
                    (self.report.loc[name]['[Adapter removal] # of reads remaining'] -
                    self.report.loc[name]['[rRNA removal] # of reads remaining']) /
                    self.report.loc[name]['[Adapter removal] # of reads remaining'] * 100)
        else:
            self.report.at[name, '[rRNA removal] # of reads remaining'] = (
                self.report.loc[name]['[Adapter removal] # of reads remaining'])
            self.report.loc[name]['[rRNA removal] % of reads removed'] = 0.0
        
        print('Quality trimming')
        # Quality trimming
        self.info_from_fastqc(output_dir, name, '[Before quality trimming]', prefix2terms)
        self.report.at[name, '[Quality trimming] Parameters'] = '; '.join(set(open(
                '{}/Preprocess/Trimmomatic/{}_quality_params.txt'.format(output_dir, 
                 name)).read().split('\n')[1:]))                                # TODO - because '' must be interfering, try to cut the problem at the root before troubles
        self.report.at[name, '[Quality trimming] # of reads remaining'] = (
            mtools.count_on_file('@', '{}/Preprocess/Trimmomatic/quality_trimmed_{}_forward_paired.fq'.format(
                    output_dir, name)))
        self.report.at[name, '[Quality trimming] % of reads removed'] = (
                (self.report.loc[name]['[rRNA removal] # of reads remaining'] -
            self.report.loc[name]['[Quality trimming] # of reads remaining']) /
            self.report.loc[name]['[rRNA removal] # of reads remaining'] * 100)
        self.info_from_fastqc(output_dir, name, '[After quality trimming]', prefix2terms)
        self.report.to_csv(output_dir + '/report.tsv', sep = '\t')
        
    def set_samples(self, sample2name):
        name2sample = {vx : k for k, v in sample2name.items() for vx in v}
        name2sample = pd.DataFrame.from_dict(name2sample, orient='index', 
                                             columns = ['Sample'])
        self.report = pd.merge(name2sample, self.report, left_index = True, 
                               right_index = True)
    
    def info_from_assembly(self, output_dir, sample):
        qc_report = pd.read_csv('{}/Assembly/{}/quality_control/report.tsv'.format(
                output_dir, sample), sep = '\t', index_col = 0).transpose()
        qc_report['Sample'] = sample
        self.report = pd.merge(self.report, qc_report, on = 'Sample', how = 'outer')
        cols = qc_report.colums.tolist(); cols.remove('Sample')
        for col in cols:
            self.report['[Assembly] ' + col] = self.report[col].fillna(
                    self.report['[Assembly] ' + col])
        self.report.drop(cols, axis=1)
        
    def info_from_annotation(self, output_dir, sample):
        sample_report = {'# of proteins detected': mtools.count_on_file('>',
                         '{}/Annotation/{}/fgs.faa'.format(output_dir, sample))}
        sample_report['# of proteins annotated (DIAMOND)'] = (sample_report['# of proteins detected'] - 
                      mtools.count_on_file('>', '{}/Annotation/{}/unaligned.fasta'.format(output_dir, sample)))
        sample_report['% of proteins annotated (DIAMOND)'] = (sample_report['# of proteins detected'] /
                      sample_report['# of proteins annotated (DIAMOND)'] * 100)
        sample_report['# of proteins annotated (PSI-BLAST)'] = mtools.count_lines(
                      '{}/Annotation/{}/cdd_aligned.blast'.format(output_dir, sample))
        sample_report['% of proteins annotated (PSI-BLAST)'] = (sample_report['# of proteins detected'] /
                      sample_report['# of proteins annotated (PSI-BLAST)'] * 100)
        tax_data = pd.read_csv()
        for col in self.taxonomic_columns:
            sample_report['Main taxa identified (' + col + ')'] = ';'.join()
            ['[Annotation] ' + col for col in (['# of proteins detected', 
         '# of proteins annotated (DIAMOND)', '% of proteins annotated (DIAMOND)',
         '# of proteins annotated (PSI-BLAST)', '% of proteins annotated (PSI-BLAST)'] +
         ['Main taxa identified (' + col + ')' for col in self.taxonomic_columns] + 
         ['Main functional categories identified'])]
        
        
        
        