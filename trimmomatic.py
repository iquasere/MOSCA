# -*- coding: utf-8 -*-
'''
Trimmomatic python API

By João Sequeira

7th March 2017
'''

import subprocess
import os
from fastqc import FastQC

class Trimmomatic:
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def update_adapters(self):
        print('Updating adapters')
        bashCommand = 'svn export https://github.com/timflutre/trimmomatic/trunk/adapters adapters'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        
    def remove_adapters(self):
        print('Determining adapters')
        input_files = [file.split('/')[-1].rstrip('.fastqc') for file in self.input_files]
        print('input_files:', input_files)
        #Check presence of adapters in input file(s)
        adapters = ['adapters/' + f for f in os.listdir('adapters') if os.path.isfile(os.path.join('adapters', f))]    
        print('adapters:', adapters)
        
        fastqc = FastQC(outdir = self.working_dir + '/preILLUMINACLIP',
                        extract = True,
                        files = self.input_files)
        fastqc.run() 
        
        adapter_contaminated = list()
        for file in input_files:
            data = self.parse_fastqc_result(self.working_dir + '/preILLUMINACLIP/' + file + '_fastqc/fastqc_data.txt')
            if not data['Overrepresented sequences'][0] == 'pass':
                no_adapters = True
                i = 0
                while i < len(data['Overrepresented sequences'][1]['Possible Source']) and no_adapters == True:
                    if 'Adapter' in data['Overrepresented sequences'][1]['Possible Source'][i]:
                        no_adapters = False
                    i += 1
                if no_adapters == False:
                    adapter_contaminated.append(file)
        
        if len(adapter_contaminated) < 0:           #no files have adapters
            return list()
        else:
            possible_adapters = list()
            for adapter in adapters:
                
                #trim according to each adapter file
                if self.paired in adapter:
                    self.illuminaclip = [adapter,'2','30','10']
                    self.output = self.working_dir + '/' +  input_files[0].split('_R')[0] + '_' + adapter.split('/')[1].rstrip('.fa')
                    self.run()

                    #generate fastqc report
                    files = [self.working_dir + '/' + input_files[0].split('_R')[0] + '_' + adapter.split('/')[1].rstrip('.fa') + '_' + fv + '_paired.fq'
                                            for fv in ['forward','reverse']]
                    fastqc = FastQC(outdir = self.working_dir,
                                    extract = True,
                                    files = files)
                    fastqc.run()
                    
                    #check presence of adapters in fastqc report                OLDES6_S4_L001_NexteraPE-PE_forward_paired.fq
                    for file in adapter_contaminated:
                        data = self.parse_fastqc_result(self.working_dir + '/' + input_files[0].split('_R')[0] + '_' + adapter.split('/')[1].rstrip('.fa') + '_forward_paired_fastqc/fastqc_data.txt')
                        if data['Overrepresented sequences'][0] == 'pass':
                            data = self.parse_fastqc_result('tool_test/trimmomatic_test/ILLUMINACLIP/' + input_files[0].split('_R')[0] + '/' + input_files[0].split('_R')[0] + '_' + adapter.split('/')[1].rstrip('.fa') + '_reverse_paired_fastqc/fastqc_data.txt')
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
                                data = self.parse_fastqc_result('tool_test/trimmomatic_test/ILLUMINACLIP/' + input_files[0].split('_R')[0] + '/' + input_files[0].split('_R')[0] + '_' + adapter.split('/')[1].rstrip('.fa') + '_reverse_paired_fastqc/fastqc_data.txt')
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
            else:
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
        
        if k == 'Per base sequence quality':            #cuts the sequences to the position before reaching red (crop)
            i = 0
            end = False
            while i < len(data['Per base sequence quality'][1]) and end == False:
                if (float(data['Per base sequence quality'][1]['Lower Quartile'][i]) < 10
                or float(data['Per base sequence quality'][1]['Median'][i]) < 25):
                    self.crop = data['Per base sequence quality'][1].index[i]
                    if '-' in self.crop:
                        self.crop = self.crop.split('-')[0]
                    end = True
                    print(self.crop)
                    self.crop = str(int(self.crop) - 1)
                i += 1
        
        if self.data != 'rrna':
            if k == 'Per base sequence content':            #cuts the sequences from the beggining until abs(A-T) and abs(G-C) are < 10
                i = 0
                while (i < len(data['Per base sequence content'][1]) and 
                    (abs(float(data['Per base sequence content'][1]['A'][i]) - float(data['Per base sequence content'][1]['T'][i])) > 10 or
                    abs(float(data['Per base sequence content'][1]['G'][i]) - float(data['Per base sequence content'][1]['C'][i])) > 10)):
                    i += 1
                self.headcrop = data['Per base sequence content'][1].index[i]
                if '-' in self.headcrop:
                    self.headcrop = self.headcrop.split('-')[1]
                self.headcrop = str(int(float(self.headcrop)))
    
    def define_by_report(self, fastqc_result):
        data = self.parse_fastqc_result(fastqc_result)
        for k,v in data.items():
            if v[0] in ['warn', 'fail']:
                self.add_fastqc_argument(data,k)
        self.run()
        
    def bash_command(self):
        result = 'java -jar ' + self.directory + 'trimmomatic.jar ' + self.paired
        if 'quality_score' in self.__dict__.keys():
            result += " -" + self.quality_score
            self.__dict__.pop('quality_score')
        for f in self.input_files:
            result += ' ' + f
        self.__dict__.pop('input_files')
        if self.paired == 'PE':
            result += ' ' + self.output + '_forward_paired.fq ' + self.output + '_forward_unpaired.fq ' + self.output + '_reverse_paired.fq ' + self.output + '_reverse_unpaired.fq'
        elif self.paired == 'SE':
             result += ' ' + self.output + '.fq'
        self.__dict__.pop('paired'); self.__dict__.pop('output')
        for arg in self.__dict__.keys():
            result += self.set_argument(arg)
        return result

    def run(self):
        bashCommand = self.bash_command()
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error