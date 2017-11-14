# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 18:33:52 2017

@author: Asus
"""

import subprocess

class Annotating:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def gene_calling(self):
        bashCommand = 'perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome=' 
        fgs_dict = {'megahit':'final.contigs.fa','metaspades':'contigs.fasta'}
        bashCommand += self.out_dir + '/Assembly/' + fgs_dict[self.assembler]
        bashCommand += ' -out=' + self.out_dir + '/Annotation/FGS -complete=1 -train=./complete'
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def gene_call(self, reads, out):
        bashCommand = 'perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome=' + reads + ' -out=' + out + ' -complete=0 -train=./illumina_10'
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
        
    def depth_of_coverage(self):
        #BEDtools
        pass
    
    def annotation(self):
        from diamond import DIAMOND
        
        diamond = DIAMOND(threads = '6',
                          db = self.db,
                          out = 'Annotation/DIAMOND/' + self.file + '/NR_MetaSPAdes',
                          query = 'Annotation/FGS.faa',
                          un = 'Annotation/Unaligned',
                          unal = '1',
                          max_target_seqs = '1')
        
        import os
        if not os.path.isfile(self.db):
            diamond.set_database(self.db, self.out_dir + '../Databases/DIAMOND')
            
        #diamond.run()
        #print(diamond.validate_annotation())
        print('lets parse')
        ids = diamond.parse_result()['sseqid']
        return ids
    
    def binning(self):
        #VizBin... duh...
        pass
    
    def run(self):
        pass