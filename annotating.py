# -*- coding: utf-8 -*-
"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

from mosca_tools import MoscaTools
mt = MoscaTools()

class Annotation:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def gene_calling(self):
        bashCommand = 'perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome='
        if self.assembled == True:
            fgs_dict = {'megahit':'final.contigs.fa','metaspades':'contigs.fasta'}
            bashCommand += self.out_dir + '/Assembly/' + fgs_dict[self.assembler]
            bashCommand += ' -out=' + self.out_dir + '/Annotation/FGS -complete=1 -train=./complete'
        else:
            fastq2fasta_command = ("paste - - - - < " + self.out_dir + "/Preprocess/SortMeRNA/rejected.fastq.fastq"
                                   + "| cut -f 1,2 | sed 's/^@/>/' | tr \"\t" "\n\" > " + self.out_dir + 
                                   "/Preprocess/SortMeRNA/rejected.fasta")
            mt.run_command(fastq2fasta_command)
            bashCommand += (self.out_dir + '/Preprocess/SortMeRNA/rejected.fasta -out=' + self.out_dir + 
            '/Annotation/FGS -complete=0 -train=./' + self.error_model)
        mt.run_command(bashCommand)
        
    def depth_of_coverage(self):
        #BEDtools
        pass
    
    def annotation(self):
        from diamond import DIAMOND
        
        diamond = DIAMOND(threads = '6',
                          db = self.db,
                          out = self.out_dir + '/Annotation/aligned.blast',
                          query = self.out_dir + '/Annotation/FGS.faa',
                          un = self.out_dir + '/Annotation/Unaligned',
                          unal = '1',
                          max_target_seqs = '1')
        
        import os
        if not os.path.isfile(self.db):
            diamond.set_database(self.db, self.out_dir + '../Databases/DIAMOND')
            
        diamond.run()
        print('lets parse')
        ids = diamond.parse_result()['sseqid']
        
        return ids
    
    def binning(self):
        #VizBin... duh...
        pass
    
    def run(self):
        self.gene_calling()
        return self.annotation()
    
def merge_fastas(file1, file2, output):
    f1 = open(file1).readlines()
    f2 = open(file2).readlines()
    out = open(output,'w')
    for i in range(0,len(f1),2):
        out.write(f1[i])
        out.write(f1[i+1])
        out.write(f2[i])
        out.write(f2[i+1])