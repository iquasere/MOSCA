# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 14:33:12 2017

@author: Asus
"""

import subprocess, os
from mosca_tools import MoscaTools

mt = MoscaTools()

class Assembling:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def set_argument(self, x):
        if isinstance(self.__dict__[x], str): 
            return ' --' + x.replace('_','-') + ' ' + self.__dict__[x]
        elif isinstance(self.__dict__[x], list): 
            result = ' --' + x.replace('_','-') + ' '
            for part in self.__dict__[x]:
                result += self.__dict__[x] + ','
            return result.rstrip(',')
        elif self.__dict__[x] == True:
            return ' --' + x.replace('_','-')
        return 'Not a valid argument'
        
    def metaspades_command(self):
        assembler = self.__dict__.pop('assembler')
        result = 'python ../../../home/jsequeira/SPAdes-3.11.1-Linux/bin/metaspades.py'
        result += ' -o ' + self.out_dir + '/Assembly'
        out_dir = self.__dict__.pop('out_dir')
        # Input data
        if hasattr(self, 'interleaved'):
            result += ' --12 ' + self.interleaved
        if hasattr(self, 'forward_paired'):
            result += ' -1 ' + self.forward_paired
        if hasattr(self, 'reverse_paired'):
            result += ' -2 ' + self.reverse_paired
        if hasattr(self, 'unpaired'):
            result += ' -s ' + self.unpaired
        data_types = {'files_paired':'-12','forward_paired':'-1','reverse_paired':'-2','unpaired':'-s','orientation':'-'}
        if hasattr(self, 'pe_libraries'):
            for library, data in self.pe_libraries.items():
                result += ' --pe' + library + data_types[data[0]]
                if data[0] != 'orientation':
                    result += ' '
                result += data[1]
        if hasattr(self, 'se_libraries'):
            for library, data in self.se_libraries.items():
                result += ' --s' + library + ' ' + data[1]
        for attr in ['interleaved', 'forward_paired', 'reverse_paired', 'unpaired']:
            if hasattr(self, attr):
                self.__dict__.pop(attr)
        for arg in self.__dict__.keys():
            result += self.set_argument(arg)
        self.out_dir = out_dir
        self.assembler = assembler
        return result
    
    def megahit_command(self):
        assembler = self.__dict__.pop('assembler')
        result = '../../../home/jsequeira/megahit/megahit -f'
        if hasattr(self, 'forward_paired') and hasattr(self, 'reverse_paired'):
            result += ' -1 ' + self.forward_paired + ' -2 ' + self.reverse_paired
            self.__dict__.pop('forward_paired'); self.__dict__.pop('reverse_paired')
        if hasattr(self, 'interleaved'):
            result += ' --12 '
            for file in self.files_paired:
                result += file + ','
            result.rstrip(',')
            self.__dict__.pop('interleaved')
            
        for arg in self.__dict__.keys():
            result += self.set_argument(arg)
        self.__dict__['assembler'] = assembler
        return result
    
    def run_assembler(self):
        if self.assembler == 'metaspades':
            bashCommand = self.metaspades_command()
        elif self.assembler == 'megahit':
            bashCommand = self.megahit_command()
        print('bash_command:', bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def run_tool(self, bashCommand):
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
    
    def bowtie2(self, reads, contigs, temp, sam, log):
        bashCommand = 'bowtie2-build ' + contigs + ' ' + temp
        #self.run_tool(bashCommand)
        bashCommand = 'bowtie2 -a -x ' + temp + ' -q -1 ' + reads[0] + ' -2 ' + reads[1]
        bashCommand = bashCommand.rstrip(',') + ' --very-sensitive -a --reorder -p 6 1> ' + sam + ' 2> ' + log
        self.run_tool(bashCommand)
        return self.parse_bowtie2(log)
        
    def parse_bowtie2(self, file):
        handler = open(file)
        line = handler.readline()
        while line:
            old = line
            line = handler.readline()
        handler.close()
        return old.split('%')[0]
    
    def metaquast(self, contigs, out_dir):
        bashCommand = 'metaquast.py --threads 6 --output-dir ' + out_dir + '/Assembly ' + contigs
        self.run_tool(bashCommand)
     
    def quality_control(self):
        from shutil import copyfile
        assembler = 'metaspades'
        terminations = {'megahit':'/final.contigs.fa', 'metaspades':'/contigs.fasta'}
        contigs = self.out_dir + '/Assembly' + terminations[self.assembler]
        temp = self.out_dir + '/Assembly/quality_control'
        sam = self.out_dir + '/Assembly/quality_control/library.sam'
        log = self.out_dir + '/Assembly/quality_control/bowtie.log'
        #self.metaquast(contigs, self.out_dir + '/Assembly/quality_control')
        percentage_of_reads = self.bowtie2([self.forward_paired,self.reverse_paired], contigs, temp, sam, log)
        
        
        if os.path.isdir(self.out_dir + '/Assembly/quality_control/combined_reference/report.tsv'):
            copyfile(self.out_dir + '/Assembly/quality_control/combined_reference/report.tsv', self.out_dir + '/quality_control/report.tsv')
        
        handler = open(self.out_dir + '/Assembly/quality_control/report.tsv', 'a')
        handler.write('Reads aligned (%)\t' + percentage_of_reads)
        
    def get_coverage_megahit(self, contigs, reads1, reads2):
        command = 'simulatedMegahit/Assembly/quality_control/library.sam'
        mt.run_command()
        
    def make_contigs_gff(self, file):
        pass
        
    def run(self):
        #self.run_assembler()
        self.quality_control()
