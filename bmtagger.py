# -*- coding: utf-8 -*-
"""
MOSCA's BMTagger wrapper for future use 
in the removal of hosts sequences

By João Sequeira

Mar 2017
"""
                                             
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