# -*- coding: utf-8 -*-
'''
DIAMOND wrapper

By João Sequeira

July 2017
'''

from mosca_tools import MoscaTools
import subprocess, os

mtools = MoscaTools()

class DIAMOND:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def set_database(self, db, output):
        bashCommand = 'diamond makedb --in ' + db + ' -d ' + output
        mtools.run_command(bashCommand)
    
    def set_argument(self, arg):
        if isinstance(self.__dict__[arg], str):
            return ' --' + arg.replace('_','-') + ' ' + self.__dict__[arg]
        return ' --' + arg.replace('_','-')
    
    def bashCommand(self):
        result = os.path.expanduser('~/diamond') + ' blastp '
        #write code for checking if .dmnd exists!!!
        for arg in self.__dict__.keys():
            result += self.set_argument(arg)
        return result
    
    def run(self):
        bashCommand = self.bashCommand()
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def validate_annotation(self):
        aligned = len(open(self.out).readlines())
        unaligned = len(open(self.un).readlines())
        return (aligned, unaligned / 2, aligned / (unaligned / 2))
    
    def parse_result(self):
        import pandas as pd
        result = pd.read_csv(self.out, sep='\t')
        result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        return result
