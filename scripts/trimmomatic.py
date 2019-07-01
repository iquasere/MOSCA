# -*- coding: utf-8 -*-
'''
Trimmomatic wrapper

By João Sequeira

March 2017
'''

from mosca_tools import MoscaTools
from fastqc import FastQC
import glob, os

mtools = MoscaTools()

class Trimmomatic:
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs

    def remove_adapters(self):
        input_files = [file.split('/')[-1].split('.fastq')[0] for file in self.input_files]
        terms_list = ['Adapter','Illumina','Primer']
        
        #Check presence of adapters in input file(s)
        adapters = (glob.glob('MOSCA/Databases/illumina_adapters/*.fa') +       # Trimmomatic's files come as .fa
                    glob.glob('MOSCA/Databases/illumina_adapters/*.fasta'))
        print('Available adapter files:', adapters)
        adapter_contaminated = list()
        
        for file in input_files:
            data = self.parse_fastqc_result(self.working_dir + '/Preprocess/FastQC/' + file + '_fastqc/fastqc_data.txt')
            no_adapters = True
            if not data['Overrepresented sequences'][0] == 'pass':
                i = 0
                while i < len(data['Overrepresented sequences'][1]['Possible Source']) and no_adapters == True:
                    for term in terms_list:
                        if term in data['Overrepresented sequences'][1]['Possible Source'][i]:
                            no_adapters = False
                    i += 1
                if no_adapters == False:
                    adapter_contaminated.append(file)
        if len(adapter_contaminated) == 0:           #no files have adapters
            print('No adapter contamination detected.')
            return []
        else:
            possible_adapters = list()
            for adapter in adapters:
                print(adapter)
                #trim according to each adapter file
                if self.paired in adapter:
                    self.illuminaclip = [adapter,'2','30','10']
                    name_n_adapter = self.name + '_' + adapter.split('/')[-1].split('.fa')[0]
                    self.output = self.working_dir + '/Preprocess/Trimmomatic/' + name_n_adapter
                    self.run()
                    #generate fastqc report
                    files = [self.working_dir + '/Preprocess/Trimmomatic/' + name_n_adapter
                                + '_' + fr + '_paired.fq' for fr in ['forward','reverse']]
                    fastqc = FastQC(outdir = self.working_dir + '/Preprocess/FastQC',
                                    extract = True,
                                    files = files)
                    fastqc.run()
                    
                    #check presence of adapters in fastqc report
                    for file in adapter_contaminated:
                        data = self.parse_fastqc_result(self.working_dir + '/Preprocess/FastQC/' + name_n_adapter  + '_forward_paired_fastqc/fastqc_data.txt')
                        if data['Overrepresented sequences'][0] == 'pass':
                            data = self.parse_fastqc_result(self.working_dir + '/Preprocess/FastQC/' + name_n_adapter  + '_reverse_paired_fastqc/fastqc_data.txt')
                            if data['Overrepresented sequences'][0] == 'pass':
                                possible_adapters.append(adapter)
                            else:
                                no_adapters = True
                                i = 0
                                while i < len(data['Overrepresented sequences'][1]['Possible Source']) and no_adapters == True:
                                    if 'Adapter' in data['Overrepresented sequences'][1]['Possible Source'][i]:
                                        no_adapters = False
                                    i += 1
                                if no_adapters == True:
                                    possible_adapters.append(adapter)
                        else:
                            no_adapters = True
                            i = 0
                            while i < len(data['Overrepresented sequences'][1]['Possible Source']) and no_adapters == True:
                                if 'Adapter' in data['Overrepresented sequences'][1]['Possible Source'][i]:
                                    no_adapters = False
                                i += 1
                            if no_adapters == True:
                                if data['Overrepresented sequences'][0] == 'pass':
                                    possible_adapters.append(adapter)
                                else:
                                    no_adapters = True
                                    i = 0
                                    while i < len(data['Overrepresented sequences'][1]['Possible Source']) and no_adapters == True:
                                        if 'Adapter' in data['Overrepresented sequences'][1]['Possible Source'][i]:
                                            no_adapters = False
                                        i += 1
                                    if no_adapters == True:
                                        possible_adapters.append(adapter)
        print('The following files might have the adapters used')
        for adapter in possible_adapters:
            print(adapter)
        return possible_adapters
        
    def parse_fastqc_result(self, file):
        import pandas as pd
        import numpy as np
        data = dict()
        fi = open(file)
        f = list()
        for line in fi:
            f.append(line.rstrip('\t\n'))
        i = 1
        while i < len(f):
            if f[i].startswith('>>') and f[i] != '>>END_MODULE':
                name, flag = f[i][2:].split('\t')[0], f[i][2:].split('\t')[1]
                if name == 'Sequence Duplication Levels':
                    i += 1
                i += 1
                labels = f[i][1:].split('\t')
                i += 1
                partial_data = np.array(labels)
                while i < len(f) and not f[i].startswith('>>'):
                    partial_data = np.append(partial_data,f[i].split('\t'))
                    i += 1
                partial_data = np.reshape(partial_data,(int(partial_data.size/len(labels)),len(labels)))
                data[name] = (flag, pd.DataFrame(data = partial_data[1:,1:],
                                                index = partial_data[1:,0],
                                                columns = partial_data[0,1:]))
            i += 1
        return data

    def set_argument(self,arg):
        if isinstance(self.__dict__[arg], str):
            return ' ' + arg.upper() + ':' + self.__dict__[arg]
        elif isinstance(self.__dict__[arg], list):
            result = ' ' + arg.upper()
            for value in self.__dict__[arg]:
                result += ":" + value
            return result
    
    def add_fastqc_argument(self,data,k):
        if k == 'Per base sequence quality':            #cuts the sequences to the position before reaching yellow (crop)
            i = 0
            end = False
            while i < len(data['Per base sequence quality'][1]) and end == False:
                if (float(data['Per base sequence quality'][1]['Lower Quartile'][i]) < 10
                or float(data['Per base sequence quality'][1]['Median'][i]) < 25):
                    crop = data['Per base sequence quality'][1].index[i]
                    if '-' in crop:
                        crop = crop.split('-')[0]
                    end = True
                    if hasattr(self, 'crop'):
                        if int(crop) < int(self.crop):
                            self.crop = crop
                            print('CROP', self.crop)
                    else:
                        self.crop = crop
                        print('CROP', self.crop)
                else:
                    i += 1
                    
        if not hasattr(self, 'data') or (hasattr(self, 'data') and self.data != 'rrna'):
            if k == 'Per base sequence content':            #cuts the sequences from the beggining until abs(A-T) and abs(G-C) are < 10
                i = int(len(data['Per base sequence content'][1]) / 2)               #don't want to cut over half of the sequences
                while i > 0 and (abs(float(data['Per base sequence content'][1]['A'][i]) - float(data['Per base sequence content'][1]['T'][i])) < 10 and
                        abs(float(data['Per base sequence content'][1]['G'][i]) - float(data['Per base sequence content'][1]['C'][i])) < 10):
                    i -= 1
                headcrop = data['Per base sequence content'][1].index[i]
                if '-' in headcrop:
                    headcrop = headcrop.split('-')[1]
                if hasattr(self, 'headcrop'):
                    if int(headcrop) > int(self.headcrop):
                        self.headcrop = str(int(headcrop) + 1)
                        print('HEADCROP', self.headcrop)
                else:
                    self.headcrop = str(int(headcrop) + 1)
                    print('HEADCROP', self.headcrop)
    
    '''
    self.input_files will be handled with arguments defined by the list of reports.
    The harshest arguments will be selected.
    '''
    def define_by_reports(self, reports):
        for report in reports:
            data = self.parse_fastqc_result(report)
            for key in ['Per base sequence quality', 'Per base sequence content']:
                if data[key][0] in ['warn', 'fail']:
                    self.add_fastqc_argument(data,key)
        self.run()
        
    def bash_command(self):
        result = ('java -jar ' + os.path.expanduser('~/anaconda3/jar/trimmomatic.jar ') 
            + self.paired + ' -threads ' + self.threads)
        if hasattr(self, 'quality_score') and self.quality_score is not None: 
            result += " -" + self.quality_score
            self.__dict__.pop('quality_score') 
        for f in self.input_files:
            result += ' ' + f
        if self.paired == 'PE':
            result += ' ' + self.output + '_forward_paired.fq ' + self.output + '_forward_unpaired.fq ' + self.output + '_reverse_paired.fq ' + self.output + '_reverse_unpaired.fq'
        elif self.paired == 'SE':
             result += ' ' + self.output + '.fq'
        # The attributes for Trimmomatic's tools work differently from the specified in problem_attributes
        problem_attributes = ['input_files', 'name', 'data', 'working_dir', 'paired', 'output', 'threads']
        attributes = dict()
        for attr in problem_attributes:
            attributes[attr] = getattr(self, attr)
        for key in problem_attributes:
            self.__dict__.pop(key)
        for arg in self.__dict__.keys():
            if self.__dict__[arg] is not None:
                result += self.set_argument(arg)
        for attr in problem_attributes:
            setattr(self, attr, attributes[attr])
        return result
    
    def run(self):
        mtools.run_command(self.bash_command())