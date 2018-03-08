# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 20:22:48 2017

@author: Asus
"""

'''
combined_reference is total
not_aligned is what suggests
summary shows for aligned organisms
'''

class MetaQUAST:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
    
    def set_argument(self,arg):
        if isinstance(self.__dict__[arg], str): 
            return ' --' + arg.replace('_','-') + ' ' + self.__dict__[arg]
        elif isinstance(self.__dict__[arg], list): 
            result = ' --' + arg.replace('_','-') + ' '
            for part in self.__dict__[arg]:
                result += self.__dict__[arg] + ','
            return result.rstrip(',')
        elif self.__dict__[arg] == True:
            return ' --' + arg.replace('_','-')
        return result
    
    def parse_report(self, file):
        data = dict()
        handler = open(file)
        handler.readline()
        line = handler.readline()
        while line:
            data[line.split('\t')[0]] = line.split('\t')[1]
            line = handler.readline().rstrip('\n')
        handler.close()
        return data
    
    def count_reads(self, file):
        lines = 0
        handler = open(file)
        for line in handler:
            lines += 1
        handler.close()
        return lines / 4
    
    
    def used_reads(self, reads, file, assembler):
        import re
        initial_reads = 0
        for reads_file in reads:
            initial_reads += self.count_reads(reads_file)
        handler = open(file)
        if assembler == 'megahit':
            line = handler.readline()
            
            m = re.search('Local assembling k = ([0-9]).+', file, flag = re.DOTALL)
    
    
    def write_report(self, files, assembler):
        
        import os
        
        parameters = ['# contigs (>= 0 bp)','# contigs (>= 1000 bp)','# contigs (>= 5000 bp)','# contigs (>= 10000 bp)',
                        '# contigs (>= 25000 bp)','# contigs (>= 50000 bp)','Total length (>= 0 bp)',
                        'Total length (>= 1000 bp)','Total length (>= 5000 bp)','Total length (>= 10000 bp)',
                        'Total length (>= 25000 bp)','Total length (>= 50000 bp)','# contigs','Largest contig',
                        'Total length','Reference length','N50','N75','L50','L75','# misassemblies',
                        '# misassembled contigs','Misassembled contigs length','# local misassemblies',
                        '# unaligned mis. contigs','# unaligned contigs','Unaligned length','Genome fraction (%)',
                        'Duplication ratio',"# N's per 100 kbp",'# mismatches per 100 kbp','# indels per 100 kbp',
                        'Largest alignment','Total aligned length']
        
        tdata = dict()
        
        for file in files:
            
            #check if has found a reference
            if os.path.isdir(os.getcwd() + '/Assembly/MetaQUAST/' + assembler + '/' + file + '/combined_reference'):
                data = self.parse_report(os.getcwd() + '/Assembly/MetaQUAST/' + assembler + '/' + file + '/combined_reference/report.tsv')
            else:
                data = self.parse_report(os.getcwd() + '/Assembly/MetaQUAST/' + assembler + '/' + file + '/report.tsv')
            for key in parameters:
                if not key in data.keys():
                    data[key] = '-'
                    
            tdata[file] = data
                    
        import pandas as pd
            
        df = pd.DataFrame.from_dict(tdata)
        
        writer = pd.ExcelWriter('assembly_quality.xlsx', engine='xlsxwriter')
        df.to_excel(writer, 'MetaSPAdes')
        writer.save()

    
    def bash_command(self):
        result = 'metaquast.py '
        file = self.input_file;self.__dict__.pop('input_file')
        if hasattr(self, 'references'):
            result += ' -R '
            for file in self.reference:
                result += file + ','
            result.rstrip(',')
            self.__dict__.pop('references')
        if hasattr(self, 'assembly_names') and self.assembly_names == True:
            result += ' -L'
            self.__dict__.pop('assembly_names')
        for arg in self.__dict__.keys():
            result += self.set_argument(arg)
        result += ' ' + file
        return result
        
    def run(self):
        import subprocess
        bashCommand = self.bash_command()
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error