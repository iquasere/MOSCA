# -*- coding: utf-8 -*-
'''
SortMeRNA python API

By João Sequeira

7th March 2017
'''

import subprocess

class SortMeRNA:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def set_optional_argument(self,arg,result):
        if hasattr(self,arg):
            result += ' --' + arg + ' ' + self.arg
        return result
    
    def generate_index(self, database):
        result = 'indexdb_rna --ref ' + database + '.fasta,' + database + '.idx'
        bashCommand = result
        process = subprocess.Popen(bashCommand.split())
        output, error = process.communicate()

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
        temp = self.working_dir + 'Preprocess/SortMeRNA/joined.fastq'
        bashCommand = 'bash MOSCA/merge-paired-reads.sh ' + self.reads[0] + ' ' + self.reads[1] + ' ' + temp
        self.reads = temp
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
    
    def divide_temp(self):
        readf, readr = self.working_dir + 'Preprocess/SortMeRNA/forward_' + self.other.split('/')[-1], self.working_dir + '/Preprocess/SortMeRNA/reverse_' + self.other.split('/')[-1]
        bashCommand = 'bash MOSCA/unmerge-paired-reads.sh ' + self.other + '.fastq ' + readf + ' ' + readr
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate() 
        
    def run_tool(self):
        bashCommand = self.bash_command()
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def run(self):
        if self.paired == True:
            self.join_fr()
            self.run_tool()
            self.divide_temp()
        else:
            self.run_tool()
 