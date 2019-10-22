# -*- coding: utf-8 -*-
"""
MOSCA's Assembly package for performing Assembly with MetaSPAdes
and Megahit and Quality Control analysis of the resulting contigs

By João Sequeira

Jun 2017
"""

import os
from mosca_tools import MoscaTools

mtools = MoscaTools()

class Assembler:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        if hasattr(self, 'memory') and self.assembler == 'metaspades':          # Megahit reads max memory in byte, metaspades reads in Gb
            self.memory /= 1.25e8
        
    def set_argument(self, x, assembler):
        if isinstance(self.__dict__[x], str):
            if x != 'threads':
                return ' --' + x.replace('_','-') + ' ' + self.__dict__[x]
            else:
                return (' --threads ' + self.__dict__[x] if assembler == 'metaspades' 
                        else ' --num-cpu-threads ' + self.__dict__[x])
        elif isinstance(self.__dict__[x], list): 
            result = ' --' + x.replace('_','-') + ' '
            for part in self.__dict__[x]:
                result += self.__dict__[x] + ','
            return result.rstrip(',')
        elif self.__dict__[x] == True:
            return ' --' + x.replace('_','-')
        else:
            print('Not a valid {} value for the assembler:{}'.format(x, self.__dict__[x]))
        
    def metaspades_command(self):
        assembler = self.__dict__.pop('assembler')
        result = 'metaspades.py'
        result += ' -o ' + self.out_dir
        out_dir = self.__dict__.pop('out_dir')
        forward, reverse = self.forward, self.reverse
        # Input data
        if hasattr(self, 'interleaved'):
            result += ' --12 ' + self.interleaved
        if hasattr(self, 'forward'):
            result += ' -1 ' + self.forward
        if hasattr(self, 'reverse'):
            result += ' -2 ' + self.reverse
        if hasattr(self, 'unpaired'):
            result += ' -s ' + self.unpaired
        data_types = {'files_paired':'-12','forward':'-1','reverse':'-2',
                      'unpaired':'-s','orientation':'-'}
        if hasattr(self, 'pe_libraries'):
            for library, data in self.pe_libraries.items():
                result += ' --pe' + library + data_types[data[0]]
                if data[0] != 'orientation':
                    result += ' '
                result += data[1]
        if hasattr(self, 'se_libraries'):
            for library, data in self.se_libraries.items():
                result += ' --s' + library + ' ' + data[1]
        for attr in ['interleaved', 'forward', 'reverse', 'unpaired']:
            if hasattr(self, attr):
                self.__dict__.pop(attr)
        for arg in self.__dict__.keys():
            if self.__dict__[arg] is not None:
                result += self.set_argument(arg, 'metaspades')
        self.out_dir = out_dir
        self.assembler = assembler
        self.forward, self.reverse = forward, reverse
        return result
    
    def megahit_command(self):
        args = dict()
        for arg in ['assembler','name']:
            args[arg] = self.__dict__.pop(arg)
        result = 'megahit -f'
        if hasattr(self, 'forward') and hasattr(self, 'reverse'):
            result += ' -1 ' + self.forward + ' -2 ' + self.reverse
            args['forward'] = self.__dict__.pop('forward')
            args['reverse'] = self.__dict__.pop('reverse')
        if hasattr(self, 'interleaved'):
            result += ' --12 '
            for file in self.files_paired:
                result += file + ','
            result.rstrip(',')
            self.__dict__.pop('interleaved')
        for arg in self.__dict__.keys():
            if self.__dict__[arg] is not None:
                result += self.set_argument(arg, 'megahit')
        for arg in args.keys():
            self.__dict__[arg] = args[arg]
        return result
    
    def run_assembler(self):
        if self.assembler == 'metaspades':
            bashCommand = self.metaspades_command()
        elif self.assembler == 'megahit':
            bashCommand = self.megahit_command()
        mtools.run_command(bashCommand)
        
    def percentage_of_reads(self, file):
        handler = open(file)
        lines = handler.readlines()
        return lines[-1].split('%')[0]
    
    def run_metaquast(self, contigs, out_dir):
        bashCommand = ('metaquast.py --threads ' + self.threads + ' --output-dir ' + 
                       out_dir + ' ' + contigs)
        mtools.run_command(bashCommand)
     
    def quality_control(self):
        contigs = self.out_dir + '/contigs.fasta'
        self.run_metaquast(contigs, self.out_dir + '/quality_control')
        mtools.perform_alignment(contigs, [self.forward, self.reverse], 
                                 self.out_dir + '/quality_control/alignment', 
                                 threads = self.threads)
        percentage_of_reads = self.percentage_of_reads(self.out_dir + '/quality_control/alignment.log')
        
        if os.path.isfile(self.out_dir + '/quality_control/combined_reference/report.tsv'):  #if metaquast finds references to the contigs, it will output results to different folders
            os.rename(self.out_dir + '/quality_control/combined_reference/report.tsv', 
                      self.out_dir + '/quality_control/report.tsv')
        
        handler = open(self.out_dir + '/quality_control/report.tsv', 'a')
        handler.write('Reads aligned (%)\t' + percentage_of_reads + '\n')
        
    def run(self):
        self.run_assembler()
        if self.assembler == 'megahit':                                         # all contigs files are outputed the metaspades way
            os.rename(self.out_dir + '/final.contigs.fa', 
                      self.out_dir + '/contigs.fasta')
        self.quality_control()                                                  # TODO - metaquast uses reads for sequence variation analysis - future work should involve determining the importance of that