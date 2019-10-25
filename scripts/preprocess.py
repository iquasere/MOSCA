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
    def rrna_removal(self):
        print('Beggining rRNA sequences removal')
        ref = (['MOSCA/Databases/rRNA_databases/' + database for database in
                    ['silva-arc-16s-id95', 'silva-arc-23s-id98', 'silva-bac-16s-id90', 
                     'silva-bac-23s-id98', 'silva-euk-18s-id95', 'silva-euk-28s-id98']])
    
        sortmerna = SortMeRNA(ref = ref,
                              reads = self.files,
                              aligned = self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_accepted',
                              output_format = ['fastx'],
                              other = self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_rejected',
                              paired = True if self.paired == 'PE' else False,
                              working_dir = self.working_dir,
                              name = self.name,
                              threads = self.threads)
        sortmerna.run()
        
        self.files = ['{}/Preprocess/SortMeRNA/{}_{}.fastq'.format(self.working_dir, self.name, fr) for fr in ['forward','reverse']]
        print('rRNA sequences removal done')
    
    # Trimmomatic - removal of low quality regions and short reads
    def quality_trimming(self):
        print('Beggining quality trimming')
        
        print('Generating quality check')
        fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                        extract = True,
                        files = self.files)
        fastqc.run()
        
        trimmomatic = Trimmomatic(input_files = self.files,
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
        
        reports = [filename.replace('SortMeRNA' if self.data == 'mrna' else     # if data is mRNA, SortMeRNA will be used
                                    'Trimmomatic', 'FastQC').split('.f')[0] + 
                   '_fastqc/fastqc_data.txt' for filename in self.files]        # .f works for SortMeRNA (.fastq) and Trimmomatic (.fq) terminations
                
        trimmomatic.define_by_reports(reports, '{}/Preprocess/Trimmomatic/{}/quality_params.txt'.format(
                self.working_dir, self.name))
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

        open('{}/Preprocess/Trimmomatic/{}/adapters.txt'.format(self.working_dir, 
             self.name), 'w').write('\n'.join(adapters))

        if len(adapters) > 0:
            adapter_part = adapters[0].split('/')[-1].rstrip('.fa')
            self.files = ['{}/Preprocess/Trimmomatic/{}_{}_{}_paired.fq'.format(
                    self.working_dir, self.name, adapter_part, fr) for fr in ['forward', 'reverse']]
        
        #self.host_sequences_removal()  
        
        if self.data == 'mrna':
            self.rrna_removal()
        
        self.quality_trimming()
        self.final_quality_check()
