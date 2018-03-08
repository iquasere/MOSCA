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
        bashCommand += ' -out=' + self.out_dir + 'Annotation/FGS -complete=1 -train=./complete'
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    '''
    def annotation(self):
        bashCommand = '../../../home/jsequeira/diamond blastp --db ' + self.db + ' --query ' + self.out_dir + 'Annotation/FGS.faa --out ' + self.out_dir + 'Annotation/aligned.blast'
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    '''
    
    def gene_call(self, reads, out):
        bashCommand = "paste - - - - < " + reads + " | cut -f 1,2 | sed 's/^@/>/' | tr \"\t" "\n\" > " + reads.rstrip('fastq') + ".fa"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        bashCommand = 'perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome=' + reads.rstrip('fastq') + '.fa -out=' + out + ' -complete=0 -train=~/FGS/train/illumina_10'
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
                          out = self.out_dir + '/Annotation',
                          query = self.out_dir + 'Annotation/FGS.faa',
                          un = self.out_dir + 'Annotation/Unaligned',
                          unal = '1',
                          max_target_seqs = '1')
        
        import os
        if not os.path.isfile(self.db):
            diamond.set_database(self.db, self.out_dir + '../Databases/DIAMOND')
            
        diamond.run()
        #print(diamond.validate_annotation())
        print('lets parse')
        ids = diamond.parse_result()['sseqid']
        return ids
    
    def binning(self):
        #VizBin... duh...
        pass
    
    def run(self):
        self.gene_calling()
        return self.annotation()