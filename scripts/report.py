# -*- coding: utf-8 -*-
'''
General report construction and export

By João Sequeira

Oct 2019
'''

from mosca_tools import MoscaTools
import pandas as pd

mtools = MoscaTools()

class Report:
    
    def __init__(self, output_dir, sample2name, mg_names, mtmp_names, 
                 metatranscriptomics = True, **kwargs):
        self.__dict__ = kwargs
        
        self.fastq_columns = ['Per base sequence quality', 'Per tile sequence quality',
                 'Per sequence quality scores', 'Per base sequence content',
                 'Per sequence GC content', 'Per base N content',
                 'Sequence Length Distribution', 'Sequence Duplication Levels', 
                 'Overrepresented sequences', 'Adapter Content']
        
        self.prefix2terms = {'[Initial quality assessment]': ('', 'R1', 'R2'),
             '[Before quality trimming]': ('after_adapter_removal_' + ),
             '[After quality trimming]': }
        
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
        
        self.build_report()
        

    def initialize_report(self):
        return pd.DataFrame(columns = (['Sample', 'Name', 'Forward/reverse'] + 
        ['[Initial quality assessment] ' + col for col in self.fastq_columns] + 
        ['[Adapter removal] adapter files', '[Adapter removal] % of reads removed',
        '[rRNA removal] % of reads removed'] + 
        ['[Before quality trimming] ' + col for col in self.fastq_columns] +
        ['[Quality trimming] ' + col for col in ['HEADCROP', 'CROP', 'AVGQUAL', 'MINLEN']] +
        ['[Quality trimming] % of reads removed'] + 
        ['[After quality trimming] ' + col for col in self.fastq_columns] + 
        ['[Assembly] ' + col for col in self.metaquast_columns] +
        ['[Annotation] ' + col for col in (['# of proteins detected', 
         '# of proteins annotated (DIAMOND)', '# of proteins annotated (PSI-BLAST)'] +
         ['Main taxa identified (' + col + ')' for col in self.taxonomic_columns] + 
         ['Main functional categories identified'])] +
        ['[Binning] ' + col for col in ['# of bins', '# of useful bins',
         'average completeness', 'average contamination']] + 
        (['[Metatranscriptomics] ' + col for col in (['# of reads aligned', 
         '# of differentially expressed proteins'] + 
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
    def info_from_fastqc(self, name, prefix, performed_rrna_removal = False):
        initial_reports = [mtools.parse_fastqc_report(
        '{}/Preprocess/FastQC/{}{}_{}_fastqc/fastqc_data.txt'.format(
                prefix2terms[prefix][0],
         self.output_dir, name, fr)) for fr in prefix2terms[prefix][1:]]
        for report in initial_reports:
            for column in self.fastq_columns:
                self.report.loc[name]['{} {}'.format(prefix, column)] = (
                        )
                        
    def build_report(self):
        # Preprocessing
        for mg_name in self.mg_names:
            self.report.append(pd.Series(name = mg_name))
            self.info_from_preprocessing(mg_name)
        if metatranscriptomics:
            for mt_name in self.mtmp_names:
                self.report.append(pandas.Series(name = mt_name))
                self.info_from_preprocessing(mt_name, performed_rrna_removal = True)
        
        # Assembly
        
        
        
        