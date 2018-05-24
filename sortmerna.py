# -*- coding: utf-8 -*-
'''
SortMeRNA wrapper

By João Sequeira

March 2017
'''

from mosca_tools import MoscaTools
from progressbar import ProgressBar

mtools = MoscaTools()

class SortMeRNA:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def set_optional_argument(self,arg,result):
        if hasattr(self,arg):
            result += ' --' + arg + ' ' + self.arg
        return result
    
    def generate_index(self, database):
        mtools.run_command('indexdb_rna --ref ' + database + '.fasta,' + database + '.idx')

    def bash_command(self):
        import os.path
        result = 'sortmerna --ref '
        for file in self.ref:
            if (os.path.isfile(file + '.idx.pos_0.dat') and os.path.isfile(file + '.idx.stats') 
            and os.path.isfile(file + '.idx.kmer_0.dat') and os.path.isfile(file + '.idx.bursttrie_0.dat')): 
                result += file + '.fasta,' + file + '.idx:'
            else:
                print('Creating index for database at ' + file)
                self.generate_index(file)
                result += file + '.fasta,' + file + '.idx:'
        result = result.rstrip(':')
        result += ' --reads ' + self.reads + ' --aligned ' + self.aligned
        for out in self.output_format:
            result += ' --' + out
        if hasattr(self,'other'):
            result += ' --other ' + self.other
        if hasattr(self,'paired_in'):
            if self.paired_in == True:
                result += ' --paired_in'
        if hasattr(self,'paired_out'):
            if self.paired_out == True:
                result += ' --paired_out' 
        print('bash command:',result)
        return result
    
    def join_fr(self):
        pbar = ProgressBar()
        handler1 = open(self.reads[0]).readlines()
        handler2 = open(self.reads[1]).readlines()
        other = self.working_dir + 'Preprocess/SortMeRNA/joined.fastq'
        temp = open(other, 'w')
        print('Merging forward ' +  self.reads[0] + ' and reverse ' + self.reads[1] + ' to interleaved ' + other)
        for i in pbar(range(0, len(handler1), 4)):
            for j in range(4):
                temp.write(handler1[i + j])
            for j in range(4):
                temp.write(handler2[i + j])
        handler1.close(); handler2.close()
    
    def divide_temp(self):
        pbar = ProgressBar()
        readf = self.working_dir + 'Preprocess/SortMeRNA/forward_' + self.other.split('/')[-1] + 'fastq'
        readr = readf.replace('forward_','reverse')
        f1 = open(readf, 'w'); f2 = open(readr, 'w')
        other = open(self.other)
        print('Splitting interleaved ' + self.other + ' to forward ' + readf + ' and reverse ' + readr)
        for i in pbar(range(0, len(self.other), 8)):
            for j in range(4):
                f1.write(other[i+j])
            for j in range(4, 8):
                f2.write(other[i+j])
            
        bashCommand = 'bash MOSCA/unmerge-paired-reads.sh ' + self.other + '.fastq ' + readf + ' ' + readr
        mtools.run_command(bashCommand)
        
    def run_tool(self):
        mtools.run_command(self.bash_command())
    
    #correct number of reads per file - if unequal number of reads from forward to reverse file, it will be corrected by separation name/1,2
    def correct_files(self):
        bashCommand = ('cat ' + self.working_dir + '/Preprocess/SortMeRNA/forward_rejected.fastq ' + 
                       self.working_dir + "/Preprocess/SortMeRNA/reverse_rejected.fastq | paste - - - - | sort | tr '\t' '\n' | seqtk dropse >" + self.working_dir + 
                       '/Preprocess/reads.fastq')
        mtools.run_command(bashCommand)
    
    def run(self):
        if self.paired == True:
            self.join_fr()
            self.run_tool()
            self.divide_temp()
            self.correct_files()
        else:
            self.run_tool()
