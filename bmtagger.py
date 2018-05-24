# -*- coding: utf-8 -*-
"""
MOSCA's BMTagger wrapper for future use 
in the removal of hosts sequences

By João Sequeira

Mar 2017
"""

#../../../../../../../../../home/sequeira/miniconda3/bin/bmtagger.sh
#human genome database extracted from http://www.ensembl.org/biomart/martview/6f3498f80caaf0f4ca4f16dd0ea0ab38

#bmtool -d bmtagger/human_exome/part0.fasta -o bmtagger/human_exome/part0.bitmask -A 0 -w 18         gives C error involved with bad allocation of memory
#srprism mkindex -i bmtagger/human_exome/part0.fasta -o bmtagger/human_exome/part0.srprism -M 7168   gives error of very large database
#makeblastdb -in bmtagger/mart_export.fa -dbtype nucl
                                             
class BMTagger:
    def __init__ (self, **kwargs):
        self.attributes = kwargs
        
    def divide_database(self, database):
        lines = list()
        for line in open(database):
            lines.append(line)
        sequences = list()
        i = 0
        print('1st done')
        while i < len(lines):
            if lines[i].startswith('>'):
                sequence = [lines[i], '']
                i += 1
                while i < len(lines) and not lines[i].startswith('>'):
                    sequence[1] += lines[i]
                    i += 1
                sequences.append(sequence)          #len for human exome: 737982
        print('2nd done')
        i = 1
        j = 0
        while i < len(sequences) - 1: #len(sequences):
            file = open('human_exome/part' + str(j) + '.fasta','w')
            print('human_exome/part' + str(j) + '.fasta')
            while i < len(sequences) - 1 and i % 100000 != 0:
                file.write(sequences[i+1][0] + sequences[i+1][1])
                i += 1
            i += 1
            j += 1
            print(j)
        print("3rd's the charm!")
        
    def generate_references(self):
        import subprocess
        bashCommand = 'bmtool -d ' + self.attributes['reference'] + '.fa -o ' + self.attributes['reference'] + '.bitmask -A 0 -w 18'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        bashCommand = 'srprism mkindex -i ' + self.attributes['reference'] + '.fa -o ' + self.attributes['reference'] + '.srprism -M 7168'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        bashCommand = 'makeblastdb -in ' + self.attributes['reference'] + '.fa -dbtype nucl'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        
    def bash_command(self):
        result = 'bmtagger.sh '
        result += self.attributes['reference'] + '.bitmask ' + self.attributes['reference'] + '.srprism --T tmp'
        if self.attributes['fasta'] == True:
            result += ' --q0'
        else:
            result += ' --q1'
        if self.attributes['paired'] == True:
            result += ' -1 ' + self.attributes['files'][0] + ' -2 ' + self.attributes['files'][1]
        else:
            result += ' -1 ' + self.attributes['files'][0]
        result += ' -­‐o ' + self.attributes['output'] + '.out'
        return result
    
    def run(self):
        import subprocess
        bashCommand = self.bash_command()
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
if __name__ == '__main__':
    bmtagger = BMTagger(reference = 'bmtagger/mart_export.txt',
                 input_file = 'real_datasets/mgm4440026.3.050.upload.fna',
                 output = 'bmtagger_output')
    bmtagger.divide_database('bmtagger/mart_export.txt')
    #bmtagger.generate_references()
    #bmtagger.run() 