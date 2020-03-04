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
        self.paired_out = self.paired
        for i in range(len(self.reads)):
            if '.gz' in self.reads[i]:
                print(self.reads[i] + ' seems to be compressed. Going to be uncrompressed.')
                mtools.run_command('gunzip ' + self.reads[i])
                self.reads[i] = self.reads[i].rstrip('.gz')
                
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
                if os.path.isfile(file):
                    print('Creating index for database at ' + file)
                    self.generate_index(file)
                    result += file + '.fasta,' + file + '.idx:'
                else:                                                           # Download rRNA databases
                    print('Downloading rRNA databases to "MOSCA/Databases/rRNA_databases"')
                    mtools.run_command('svn export https://github.com/biocore/sortmerna/trunk/data/rRNA_databases MOSCA/Databases/rRNA_databases')
                    mtools.run_pipe_command('find MOSCA/Databases/rRNA_databases/* | grep -v ".fasta" | xargs rm -fr')
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
    
    '''
    Input:
        filename: str - filename of FastQ file with irregular reads (less than
                        4 lines per read)
    Output:
        Irregular reads will be removed
    '''  
    def remove_messed_reads(self, filename):
        mtools.run_pipe_command("awk 'BEGIN {RS=\"@\"; FS=\"\\n\"} { if (NF == 5) " +
        "print \"@\" substr($0,1,length-1) }' " + filename + " > " + filename + ".temp")
        os.rename(filename + ".temp", filename)
    
    # correct number of reads per file - if unequal number of reads from forward to reverse file, it will be corrected by separation name/1,2
    # from www.biostars.org/p/6925/#6928
    def remove_orphans(self, forward, reverse):
        mtools.run_pipe_command("awk '{printf substr($0,1,length-2);getline;printf \"\\t\"$0;getline;getline;print \"\\t\"$0}' "
                                       + forward + "| sort -T. > " + self.working_dir + '/Preprocess/SortMeRNA/read1.txt')
        mtools.run_pipe_command("awk '{printf substr($0,1,length-2);getline;printf \"\\t\"$0;getline;getline;print \"\\t\"$0}' "
                                       + reverse + "| sort -T. > " + self.working_dir + '/Preprocess/SortMeRNA/read2.txt')
        mtools.run_pipe_command("join " + ' '.join(["{}/Preprocess/SortMeRNA/{}".format(
                self.working_dir, fr) for fr in ['read1.txt','read2.txt']]) + 
                " | awk '{{print $1\" \"$2\"\\n\"$3\"\\n+\\n\"$4 > \"{}\";print $1\" \"$5\"\\n\"$6\"\\n+\\n\"$7 > \"{}\"}}'".format(
                                   forward, reverse))
        for file in ["{}/Preprocess/SortMeRNA/read{}.txt".format(self.working_dir, 
                     number) for number in ['1','2']]:
            os.remove(file)
    
    def run(self):
        if self.paired == True:
            basename = self.working_dir + '/Preprocess/SortMeRNA/' + self.name
            interleaved = basename + '_interleaved.fastq'
            for i in range(len(self.reads)):
                if self.reads[i].endswith('.gz'):
                    mtools.run_command('gunzip ' + self.reads[i])
                    self.reads[i] = self.reads[i].rstrip('.gz')
                    
            self.merge_pe(self.reads[0], self.reads[1], interleaved)
            self.reads = interleaved
            self.run_tool()
            self.unmerge_pe(basename + '_rejected.fastq', 
                            basename + '_forward.fastq', 
                            basename + '_reverse.fastq')
            
            for fr in ['forward', 'reverse']:
                self.remove_messed_reads('{}_{}.fastq'.format(basename, fr))
            
            self.remove_orphans(basename + '_forward.fastq', 
                               basename + '_reverse.fastq')
        else:
            self.run_tool()
