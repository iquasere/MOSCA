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

class Preprocesser:

    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
    
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
        trimmomatic = Trimmomatic(input_files = self.files,
                                  paired = self.paired,
                                  minlen = '100',
                                  working_dir = self.working_dir,
                                  output = self.working_dir + '/Preprocess/Trimmomatic/after_adapter_removal_' + self.name,
                                  data = self.data,
                                  name = self.name,
                                  threads = self.threads)
        if hasattr(self, 'quality_score'):
            setattr(trimmomatic, 'quality_score', self.quality_score)
        adapters = trimmomatic.remove_adapters()
        return adapters
        print('Adapter removal done')
            
    #SortMeRNA - taxonomic assignment
    def rrna_removal(self, reads):
        print('Beggining rRNA sequences removal')
        ref = (['MOSCA/Databases/rRNA_databases/' + database for database in
                    ['silva-arc-16s-id95', 'silva-arc-23s-id98', 'silva-bac-16s-id90', 
                     'silva-bac-23s-id98', 'silva-euk-18s-id95', 'silva-euk-28s-id98']])
        sortmerna = SortMeRNA(ref = ref,
                              reads = reads,
                              aligned = self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_accepted',
                              output_format = ['fastx'],
                              other = self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_rejected',
                              paired = True if self.paired == 'PE' else False,
                              working_dir = self.working_dir,
                              name = self.name)
        sortmerna.run()
        print('rRNA sequences removal done')
    
    def quality_trimming(self):
        print('Beggining quality trimming')
        
        print('Generating quality check')
        fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                        extract = True,
                        files = [self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_' + fr + '.fastq' for fr in ['forward','reverse']])
        fastqc.run()
        
        trimmomatic = Trimmomatic(input_files = [self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_' + fr + '.fastq' for fr in ['forward','reverse']],
                                  paired = self.paired,
                                  working_dir = self.working_dir,
                                  avgqual = '20',
                                  data = self.data,
                                  name = self.name,
                                  output = self.working_dir + '/Preprocess/Trimmomatic/quality_trimmed_' + self.name,
                                  minlen = '100',
                                  threads = self.threads)
        if hasattr(self, 'quality_score'):
            setattr(trimmomatic, 'quality_score', self.quality_score)
        reports = [self.working_dir + '/Preprocess/FastQC/' + self.name + '_' + fr + '_fastqc/fastqc_data.txt' for fr in ['forward','reverse']]
        trimmomatic.define_by_reports(reports)
        print('Quality trimming done')
        
    #BMTagger - removal of human sequences 
    # TODO - implement
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
    
    def final_quality_check(self):
        print('Beggining third quality check')
        fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                        extract = True,
                        files = [self.working_dir + '/Preprocess/Trimmomatic/quality_trimmed_' + self.name + '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']])
        fastqc.run()
        print('Third quality check done')
        
        
    def run(self):
        
        self.first_check()
        
        adapters = self.trim_adapters()
        
        print(adapters)
        
        if len(adapters) > 0:
            adapter_part = adapters[0].split('/')[-1].rstrip('.fa')
            reads = [self.working_dir + '/Preprocess/Trimmomatic/' + self.name + '_' + adapter_part + '_' + fr + '_paired.fq' for fr in ['forward', 'reverse']]
        else:
            reads = self.files
        
        self.rrna_removal(reads)
        
        self.quality_trimming()
        self.final_quality_check()
        #self.host_sequences_removal()
        
if __name__ == '__main__':
    
    for files in [['Datasets/4478-DNA-S1611-MiSeqKapa/4478-R1-1-MiSeqKapa_R1.fastq', 'Datasets/4478-DNA-S1611-MiSeqKapa/4478-R1-1-MiSeqKapa_R2.fastq']]:
        preprocesser = Preprocesser(files = files,
                                    working_dir = 'MOSCAfinal',
                                    name = files[0].split('/')[-1].split('_R')[0],
                                    paired = 'PE',
                                    trimmomatic_dir = os.path.expanduser('~/anaconda3/bin'),
                                    data = 'mrna')
        preprocesser.run()
        
