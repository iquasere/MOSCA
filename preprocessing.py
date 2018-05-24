# -*- coding: utf-8 -*-
'''
MOSCA Preprocessing package for quality check and removal of
undesired reads by quality, detection of artificial origin
or detection as rRNA

By João Sequeira

March 2017
'''

from fastqc import FastQC
from trimmomatic import Trimmomatic
from sortmerna import SortMeRNA
from bmtagger import BMTagger

import os

class Preprocessing:

    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        self.name = self.files[0].split('/')[-1].split('_R')[0]
    
    #FASTQC - quality check
    def first_check(self):
        print('Beggining first quality check')
        fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                        extract = True,
                        files = self.files)
        fastqc.run()
        print('First quality check done')

    #Trimmomatic - removal of adapters
    def trim_adapters(self):
        print('Checking for adapters')
        trimmomatic = Trimmomatic(directory = '../../../home/jsequeira/anaconda3/bin/',
                                  input_files = self.files,
                                  paired = self.paired,
                                  minlen = '100',
                                  working_dir = self.working_dir,
                                  output = self.working_dir + '/Preprocess/Trimmomatic/after_adapter_removal_' + self.name,
                                  data = self.data,
                                  name = self.name)
        adapters = trimmomatic.remove_adapters()
        return adapters
        print('Adapter removal done')
    
    def quality_trimming(self, adapters = []):
        print('Beggining quality trimming')
        print('Adapters', adapters)
        if len(adapters) > 0:
            for adapter in adapters:
                trimmomatic = Trimmomatic(directory = '../../../home/jsequeira/anaconda3/bin/',
                                          input_files = [self.working_dir + '/Preprocess/Trimmomatic/' + adapter.split('/')[-1].rstrip('.fa') + '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']],
                                          paired = self.paired,
                                          working_dir = self.working_dir,
                                          output = self.working_dir + '/Preprocess/quality_trimmed_' + self.name.rstrip('/'),
                                          avgqual = '20',
                                          data = self.data,
                                          name = self.name)
                trimmomatic.define_by_report(adapter)
        else:
            trimmomatic = Trimmomatic(directory = '../../../home/jsequeira/anaconda3/bin/',
                                      input_files = self.files,
                                      paired = self.paired,
                                      working_dir = self.working_dir,
                                      output = self.working_dir + '/Preprocess/quality_trimmed_' + self.name.rstrip('.fastq'),
                                      avgqual = '20',
                                      data = self.data,
                                      name = self.name)
            trimmomatic.define_by_report()
        print('Quality trimming done') 
        
    #BMTagger - removal of human sequences
    def host_sequences_removal(self):        
        print('Beggining host sequences removal')
        bmtagger = BMTagger(script_dir = '~/miniconda3/bin',
                            reference = 'bmtagger/mart_export.txt',
                            files = 'real_datasets/mgm4440026.3.050.upload.fna',
                            output = 'bmtagger_output',
                            fasta = True,
                            paired = 'PE')
        bmtagger.run()
        print('Host sequences removal done')
    
    #SortMeRNA - taxonomic assignment
    def rrna_removal(self):
        print('Beggining rRNA sequences removal')
        ref = (['Databases/rRNA_databases/SILVA/' + database for database in
                    ['silva-arc-16s-id95', 'silva-arc-23s-id98', 'silva-bac-16s-id90', 
                     'silva-bac-23s-id98', 'silva-euk-18s-id95', 'silva-euk-28s-id98']])
        sortmerna = SortMeRNA(ref = ref,
                              reads = [self.working_dir + '/Preprocess/quality_trimmed_' + self.name + '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']],
                              aligned = self.working_dir + '/Preprocess/SortMeRNA/accepted',
                              output_format = ['fastx'],
                              other = self.working_dir + '/Preprocess/SortMeRNA/rejected',
                              paired = True if self.paired == 'PE' else False,
                              working_dir = self.working_dir)
        sortmerna.run()
        print('rRNA sequences removal done')
    
    def final_quality_check(self):
        print('Beggining third quality check')
        fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                        extract = True,
                        files = [self.working_dir + '/Preprocess/SortMeRNA/rejected_' + self.name + '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']])
        fastqc.run()
        print('Third quality check done')
        
    def run(self):
        print(self.name)
        #self.first_check()
        #adapters = self.trim_adapters()
        adapters = ['adapters/TruSeq2-PE.fa','adapters/TruSeq3-PE-2.fa']
        print(adapters)
        #self.quality_trimming(adapters)
        self.rrna_removal()
        #self.host_sequences_removal()
