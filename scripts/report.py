# -*- coding: utf-8 -*-
'''
General report construction and export

By João Sequeira

Oct 2019
'''

from mosca_tools import MoscaTools
import pandas as pd, numpy as np, os

mtools = MoscaTools()

class Report:
    
    def __init__(self, output_dir, sample2name, mg_names, mtmp_names, r1_file,
                 r2_file = None, metatranscriptomics = True, adapter = None, **kwargs):
        self.__dict__ = kwargs
        
        self.fastq_columns = ['Per base sequence quality', 'Per tile sequence quality',
                 'Per sequence quality scores', 'Per base sequence content',
                 'Per sequence GC content', 'Per base N content',
                 'Sequence Length Distribution', 'Sequence Duplication Levels', 
                 'Overrepresented sequences', 'Adapter Content']
        
        self.prefix2terms = {'[Initial quality assessment]': ('', 'R1', 'R2'),
             '[Before quality trimming]': (('', 'forward_paired', 'reverse_paired') if metatranscriptomics 
             else ('after_adapter_removal_', '_' + adapter + '_forward_paired', 
                   '_' + adapter + '_reverse_paired') if adapter is not None else ('', 'R1', 'R2')),
             '[After quality trimming]': ('quality_trimmed_', 'forward_paired', 'reverse_paired')}
        
        self.metaquast_columns = ['# contigs (>= 0 bp)', '# contigs (>= 1000 bp)', 
                     '# contigs (>= 5000 bp)', '# contigs (>= 10000 bp)',
                     '# contigs (>= 25000 bp)', '# contigs (>= 50000 bp)',
                     'Total length (>= 0 bp)', 'Total length (>= 1000 bp)',
                     'Total length (>= 5000 bp)', 'Total length (>= 10000 bp)',
                     'Total length (>= 25000 bp)', 'Total length (>= 50000 bp)',
                     '# contigs', 'Largest contig', 'Total length', 
                     'Reference length', 'N50', 'N75', 'L50', 'L75',
                     '# misassemblies', '# misassembled contigs', 
                     'Misassembled contigs', '# local misassemblies',
                     '# unaligned mis. contigs', '# unaligned contigs',
                     'Unaligned length', 'Genome fraction (%)',
                     'Duplication ratio', "# N's per 100 kbp",
                     '# mismatches per 100 kbp', '# indels per 100 kbp',
                     'Largest alignment', 'Total aligned length', 'Reads aligned (%)']
        
        self.taxonomic_columns = ['superkingdom', 'phylum', 'class', 'order', 'family', 
                     'genus', 'species']
        
        self.name2sample = {vx : k for k, v in sample2name.items() for vx in v}
        
        self.build_report()
        

    def initialize_report(self):
        return pd.DataFrame(columns = (['Sample', '# of initial reads'] +
        ['[Initial quality assessment] ' + col for col in self.fastq_columns] + 
        ['[Adapter removal] adapter files', '[Adapter removal] # of reads remaining', 
        'Adapter removal] % of reads removed', '[rRNA removal] # of reads remaining',
        '[rRNA removal] % of reads removed'] + 
        ['[Before quality trimming] ' + col for col in self.fastq_columns] +
        ['[Quality trimming] ' + col for col in ['HEADCROP', 'CROP', 'AVGQUAL', 'MINLEN']] +
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
         ['Main functional categories identified'])] if self.metatranscriptomics else 
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
    def info_from_fastqc(self, output_dir, name, prefix):
        reports = [mtools.parse_fastqc_report(
                '{}/Preprocess/FastQC/{}{}_{}_fastqc/fastqc_data.txt'.format(output_dir, 
                 self.prefix2terms[prefix][0], name, self.prefix2terms[prefix][i]) for i in [1, 2])]
        for column in self.fastq_columns:
            if reports[0][column][0] == reports[1][column][0]:
                self.report.loc[name]['{} {}'.format(prefix, column)] = reports[0][column][0]
            else:
                self.report.loc[name]['{} {}'.format(prefix, column)] = (
                        '{} (forward) {} (reverse)'.format(
                                reports[0][column][0], reports[1][column][0]))
                        
    def info_from_preprocessing(self, output_dir, name, performed_rrna_removal = False):
        self.report.loc[name]['Sample'] = self.name2sample[name]
        self.report.loc[name]['# of initial reads'] = mtools.count_on_file(self.r1_file, '@')
        self.info_from_fastqc(output_dir, name, '[Initial quality assessment]')
        
        if os.isfile('{}/Preprocess/Trimmomatic/{}/adapters.txt'.format(output_dir, name)):
            adapter_files = [file.split('/').rstrip('.fa') for file in 
                open('{}/Preprocess/Trimmomatic/{}/adapters.txt').read().split('\n')]
            self.report.loc[name]['[Adapter removal] adapter files'] = ', '.join(adapter_files)
            self.report.loc[name]['[Adapter removal] # of reads remaining'] = (
                mtools.count_on_file('{}/Preprocess/Trimmomatic/after_adapter_removal_{}_{}_forward_paired.fq'.format(
                    output_dir, name, adapter_files[0]), '@'))
            self.report.loc[name]['[Adapter removal] % of reads removed'] = (
                    self.report.loc[name]['[Adapter removal] # of reads remaining'] /
                    self.report.loc[name]['# of initial reads'] * 100)
        else:
            self.report.loc[name]['[Adapter removal] adapter files'] = np.nan
            self.report.loc[name]['[Adapter removal] # of reads remaining'] = (
                self.report.loc[name]['# of initial reads'])
            self.report.loc[name]['[Adapter removal] % of reads removed'] = 0.0
        
        if performed_rrna_removal:
            self.report.loc[name]['[rRNA removal] # of reads remaining'] = (
            mtools.count_on_file('{}/Preprocess/SortMeRNA/{}_forward_paired.fq'.format(
                    output_dir, name), '@'))
            self.report.loc[name]['[rRNA removal] % of reads removed'] = (
                    self.report.loc[name]['[rRNA removal] # of reads remaining'] /
                    self.report.loc[name]['[Adapter removal] # of reads remaining'] * 100)
        else:
            self.report.loc[name]['[rRNA removal] # of reads remaining'] = (
                self.report.loc[name]['[Adapter removal] # of reads remaining'])
            self.report.loc[name]['[rRNA removal] % of reads removed'] = 0.0
            
        self.info_from_fastqc(output_dir, name, '[Before quality trimming]')
        qual_params = ['HEADCROP', 'CROP', 'AVGQUAL', 'MINLEN']
        qual_params_vals = open('{}/Preprocess/Trimmomatic/{}/quality_params.txt'.format(
                output_dir, name)).read().split('\n')
        for i in range(4):
            self.report.loc[name]['[Quality trimming] ' + qual_params[i]] = qual_params_vals[i]
        self.report.loc[name]['[Quality trimming] # of reads remaining'] = (
            mtools.count_on_file('{}/Preprocess/Trimmomatic/quality_trimmed_{}_forward_paired.fq'.format(
                    output_dir, name), '@'))
        self.report.loc[name]['[Quality trimming] % of reads removed'] = (
            self.report.loc[name]['[Quality trimming] # of reads remaining'] /
            self.report.loc[name]['[rRNA removal] # of reads remaining'])
        self.info_from_fastqc(output_dir, name, '[After quality trimming]')
                        
    def build_report(self):
        # Preprocessing
        for mg_name in self.mg_names:
            self.report.append(pd.Series(name = mg_name))
            self.info_from_preprocessing(mg_name)
        if metatranscriptomics:
            for mt_name in self.mtmp_names:
                self.report.append(pd.Series(name = mt_name))
                self.info_from_preprocessing(mt_name, performed_rrna_removal = True)
        
        # Assembly
        
        
        
        