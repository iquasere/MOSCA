# -*- coding: utf-8 -*-
'''
SortMeRNA wrapper

By João Sequeira

March 2017
'''

from mosca_tools import MoscaTools

import os

mtools = MoscaTools()

class SortMeRNA:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        '''
        self.paired_out = self.paired
        for i in range(len(self.reads)):
            if '.gz' in self.reads[i]:
                print(self.reads[i] + ' seems to be compressed. Going to be uncrompressed.')
                mtools.run_command('gunzip ' + self.reads[i])
                self.reads[i] = self.reads[i].rstrip('.gz')
        '''
                
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
        result += ' -a ' + self.threads
        return result
    
    def merge_pe(self, forward, reverse, interleaved):
        mtools.run_command('bash MOSCA/scripts/merge-paired-reads.sh ' + forward + 
                           ' ' + reverse + ' ' + interleaved)
        
    def unmerge_pe(self, interleaved, forward, reverse):
        mtools.run_command('bash MOSCA/scripts/unmerge-paired-reads.sh ' + 
                           interleaved + ' ' + forward + ' ' + reverse)
        
    def run_tool(self):
        mtools.run_command(self.bash_command())
    
    # correct number of reads per file - if unequal number of reads from forward to reverse file, it will be corrected by separation name/1,2
    # from www.biostars.org/p/6925/#6928
    def correct_files(self, forward, reverse):
        mtools.run_pipe_command("awk '{printf substr($0,1,length-2);getline;printf \"\\t\"$0;getline;getline;print \"\\t\"$0}' "
                                       + forward + "| sort -S 8G -T. > " + self.working_dir + '/Preprocess/SortMeRNA/read1.txt')
        mtools.run_pipe_command("awk '{printf substr($0,1,length-2);getline;printf \"\\t\"$0;getline;getline;print \"\\t\"$0}' "
                                       + reverse + "| sort -S 8G -T. > " + self.working_dir + '/Preprocess/SortMeRNA/read2.txt')
        mtools.run_pipe_command("join " + ' '.join(["{}/Preprocess/SortMeRNA/{}".format(
                self.working_dir, fr) for fr in ['read1.txt','read2.txt']]) + 
                " | awk '{{print $1\" \"$2\"\\n\"$3\"\\n+\\n\"$4  > \"{}\";print $1\" \"$5\"\\n\"$6\"\\n+\\n\"$7 > \"{}\"}}'".format(
                                   forward, reverse))
        for file in ["{}/Preprocess/SortMeRNA/read{}.txt".format(self.working_dir, 
                     number) for number in ['1','2']]:
            os.remove(file)
    
    def run(self):
        if self.paired == True:
            interleaved = self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_interleaved.fastq'
            self.merge_pe(self.reads[0], self.reads[1], interleaved)
            self.reads = interleaved
            self.run_tool()
            self.unmerge_pe(self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_rejected.fastq', 
                            self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_forward.fastq', 
                            self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_reverse.fastq')
            self.correct_files(self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_forward.fastq', 
                               self.working_dir + '/Preprocess/SortMeRNA/' + self.name + '_reverse.fastq')
        else:
            self.run_tool()
            
if __name__ == '__main__':
    smr = SortMeRNA(working_dir = 'MOSCAfinal')
    
    smr.correct_files('MOSCAfinal/Preprocess/SortMeRNA/4478-R1-1-MiSeqKapa_forward.fastq',
                      'MOSCAfinal/Preprocess/SortMeRNA/4478-R1-1-MiSeqKapa_reverse.fastq')
    
    
