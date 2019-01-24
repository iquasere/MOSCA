# -*- coding: utf-8 -*-
"""
RNA-Seq simulator for differential expression

Created on Sun Jan 14 19:14:01 2018
"""

from annotation import Annotater
from diamond import DIAMOND
import pandas as pd
from mosca_tools import MoscaTools
from progressbar import ProgressBar
import os

annotater = Annotater()
mtools = MoscaTools()

class RNASeqSim:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        
    def use_fgs(self, genome, output):
        command = ('perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome='
                   + genome + ' -out=' + output + ' -complete=1 -train=./complete')
        mtools.run_command(command)
        
    def use_diamond(self, orfs, output):
        diamond = DIAMOND(threads = '6',
                          out = output,
                          query = orfs,
                          max_target_seqs = '1',
                          db = 'Databases/Annotation/uniprot.dmnd')
        diamond.run()
        
    def uniprot_mapping(self, ids, max_retry = 3):
        import requests, time
        
        BASE_URL = 'http://www.uniprot.org/uniprot/'
        
        params = {
            'from':'ACC+ID',
            'to':'ACC+ID',
            'format':'tab',
            'query':'+OR+'.join(['accession:'+acc for acc in ids]),
            'columns':'id,pathway,protein names'
            }
        response = requests.post(BASE_URL, params=params)
        cnt=0
        max_retry = max_retry
        while cnt < max_retry and response.status_code != requests.codes.ok:
            cnt += 1
            time.sleep(60)
            response = requests.post(BASE_URL, params=params)
        return [match.split('\t') for match in response.content.decode('utf-8').split('\n')[1:-1]]
        
    def chunky(self, ids, step = 100, max_retry = 3):
        result = list()
        for i in range(0, len(ids), step):
            query = ids[i: i + step]
            result += self.uniprot_mapping(query, max_retry)
        return pd.DataFrame([row[1:3] for row in result], index = [row[0] for row in result], columns = ['pathway','name'])
        
    def pathway_abundances(self, file):
        handler = open(file)
        lines = [line.rstrip('\n') for line in handler]
        return {line.split('\t')[0]:line.split('\t')[1] for line in lines}
    
    def parse_fgs(self, file):
        handler = open(file)
        result = pd.DataFrame(columns = ['sequence'])
        lines = [line.rstrip('\n') for line in handler]
        for line in lines:
            if line.startswith('>'):
                name = line[1:]
                result.loc[name] = ''
            else:
                result.loc[name] += line
        return result
        
    def join_information(self, genome, output):
        annotater.use_fgs(genome, output + '/fgs')
        annotater.use_diamond(output + '/fgs.faa', output + '/aligned.blast')
        blast = DIAMOND(out = output + '/aligned.blast').parse_result()
        blast = blast.drop_duplicates(subset = 'sseqid')
        ids = list(set([ide.split('|')[1] if ide != '*' else ide for ide in blast['sseqid']]))
        pathinfo = self.chunky(ids)
        blast.index = [ide.split('|')[1] if ide != '*' else ide for ide in blast['sseqid']]
        result = pd.concat([blast, pathinfo], axis=1, join='inner')
        result.index = result.qseqid
        result.to_csv(output + '/pathwayinfo', sep = '\t')
        return result[['qseqid', 'pathway','name']]

    #abundance is tuple (organism_presence, {pathway:expression})
    def define_abundance(self, excel, profiles, uniprotinfo, output, 
                         factor = 1, base = 0):
        simulated = pd.read_excel(excel, index = False)
        simulated = simulated[simulated['Transcriptome link'].notnull()].reset_index()
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t').drop_duplicates()
        uniprotinfo.columns = ['EMBL'] + uniprotinfo.columns.tolist()[1:]
        handler = open(output,'w')
        print('Defining expressions for ' + str(len(simulated)) + ' organisms.')
        for i in range(len(simulated)):
            pbar = ProgressBar()
            print('Dealing with ' + simulated.iloc[i]['Species'])
            rna_file = 'SimulatedMGMT/Transcriptomes/' + simulated.iloc[i]['Transcriptome link'].split('/')[-1].split('.gz')[0]
            rnas = mtools.parse_fasta(rna_file).keys()
            abundance = simulated.iloc[i]['Abundance']
            profile = profiles[simulated.iloc[i]['Profile']]
            for rna in pbar(rnas):
                found = False
                ide = rna.split()[0]
                info = uniprotinfo[uniprotinfo['EMBL'] == ide]
                for path in profile['Pathway'].keys():
                    if path in str(info['Pathway']): 
                        handler.write(rna + '\t' + str(base + factor * abundance) + '\n')
                        found = True
                if found == False:
                    if 'ATP synthase' in str(info['Protein names']):
                        handler.write(rna + '\t' + str(base + factor * abundance) + '\n')
                    else:
                        handler.write(rna + '\t' + str(base + abundance) + '\n')
        
    
        
    def use_grinder(self, fasta, abundance, total_reads, output):
        command = (os.path.expanduser('~/biogrinder/script/grinder') + 
                   ' -reference_file ' + fasta + ' -abundance_file ' + abundance + 
                   ' -total_reads ' + str(total_reads) + ' -mate_orientation FR -random_seed ' + 
                   '13 -fastq_output 1 -qual_levels 30 10 -output_dir ' + output + 
                   ' -read_dist 151 -insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8')
        mtools.run_command(command)
        
if __name__ == '__main__':
    import glob, pathlib, os
    
    files = glob.glob('SimulatedMGMT/Transcriptomes/*.fa')
    
    rnaseqsimer = RNASeqSim()
    '''
    annotater = Annotater()
    
    simulated = pd.read_excel('SimulatedMGMT/simulated_taxa.xlsx', index=False)
    simulated = simulated[simulated['Transcriptome link'].notnull()]
    simulated = simulated.reset_index()
    uniprotinfo = pd.read_csv('SimulatedMGMT/Transcriptomes/uniprot.info',sep='\t')
    
    abundances = {'archaea_co2': {
                        'Pathway':{
                                'One-carbon metabolism; methanogenesis from CO(2)':1, 
                                'Cofactor biosynthesis':1}, 
                        'Protein names':{
                                'ATP synthase':1}
                        },
                'archaea_acetate': {
                        'Pathway':{
                                'One-carbon metabolism; methanogenesis from acetate':1, 
                                'Cofactor biosynthesis':1}, 
                        'Protein names':{
                                'ATP synthase':1}
                        },
                'bacteria':{
                        'Pathway':{
                                'lipid metabolism': 1, 
                                'Cofactor biosynthesis':1,
                                'Metabolic intermediate biosynthesis': 1}, 
                        'Protein names':{
                                'ATP synthase':1}
                        }
                    }
                    
    for factor in [3,0.17,1]:
        rnaseqsimer.define_abundance('SimulatedMGMT/simulated_taxa.xlsx',
                           abundances, 'SimulatedMGMT/Transcriptomes/uniprot.info',
                           'SimulatedMGMT/Transcriptomes/abundance' + str(factor) + '.config',
                           factor = factor)
    '''
    
    for factor in [1,3,0.17]:
        rnaseqsimer.use_grinder('SimulatedMGMT/Transcriptomes/all_transcriptome.fa',
                              'SimulatedMGMT/Transcriptomes/abundance' + str(factor) + '.config',
                              1000000, 'SimulatedMGMT/Transcriptomes/b/' + str(factor))
    
    for factor in [1,3,0.17]:
        rnaseqsimer.use_grinder('SimulatedMGMT/Transcriptomes/all_transcriptome.fa',
                              'SimulatedMGMT/Transcriptomes/abundance' + str(factor) + '.config',
                              1000000, 'SimulatedMGMT/Transcriptomes/c/' + str(factor))
    
    
    
    '''
    simulated = pd.read_excel('SimulatedMGMT/Genomes/simulated_taxa.xlsx', index = False, header = None)
    simulated = simulated[simulated[8].str.contains('ftp')]
    simulated = simulated.reset_index()
    
    abundance_handler = open('SimulatedMGMT/Genomes/abundance.config','w')
    all_genomes_handler = open('SimulatedMGMT/Genomes/genomes.fasta','w')
    
    for i in range(len(simulated)):
        file = 'SimulatedMGMT/Genomes/' + simulated.loc[i][6].replace(' ','_') + '.fasta'
        fasta = mtools.parse_fasta(file)
        for k,v in fasta.items():
            abundance_handler.write(simulated.loc[i][6] + '\t' + str(simulated.loc[i][7]) + '\n')
            all_genomes_handler.write('>' + simulated.loc[i][6] + '\n' + v + '\n')
    
    command = ('/home/jsequeira/biogrinder/script/grinder -reference_file SimulatedMGMT/Genomes/genomes.fasta -abundance_file SimulatedMGMT/Genomes/abundance.config ' + 
                   '-total_reads 12834676 -mate_orientation FR -random_seed 13 -fastq_output 1 ' + 
                   '-qual_levels 30 10 -output_dir SimulatedMGMT/Genomes/MG_reads -read_dist 151 -insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8')
    mtools.run_command(command)
    '''
    '''
    for number in ['1']:
        abundancedf = pd.read_csv('../PostThesis/piptest/abundance' + number + '.config', sep = '\t', header = None, index_col = 0)
        simulator = RNASeqSim()
        
        pbar = ProgressBar()
        already_done = next(os.walk('../PostThesis/piptest/grinder/MT'))[1]
        files1 = glob.glob('../PostThesis/piptest/grinder/MT/*.fasta')
        files = [file for file in pbar(files1) if file.split('/')[-1].rstrip('.fasta') not in already_done]
        print(str(len(files1)) + ' ' + str(len(files)))
        
        for file in files:
            name = file.rstrip('.fasta').split('/')[-1]        
            with open('../PostThesis/piptest/grinder/MT/abundance.config','w') as f:
                f.write(name + '\t1')
            directory = '../PostThesis/piptest/grinder/MT/' + name
            path = pathlib.Path(directory)
            path.mkdir(parents=True, exist_ok=True)
            simulator.use_grinder(directory, float(abundancedf.loc[name]), name, factor = 1)
    '''
    '''    
    1223280
    pbar = ProgressBar()
    
    from mosca_tools import MoscaTools
    
    mt = MoscaTools()
    file = mt.parse_fasta('transcripts.fasta')
    
    for k in pbar(file.keys()):
        handler = open('../PostThesis/piptest/grinder/MT/' + k + '.fasta', 'w')
        handler.write('>' + k + '\n' + file[k] + '\n')
    
    archaea_co2 = {'One-carbon metabolism; methanogenesis from CO(2)':1, 
                   'Cofactor biosynthesis':1, 'ATP synthase':1}
    archaea_acetate = {'One-carbon metabolism; methanogenesis from acetate':1, 
                       'Cofactor biosynthesis':1, 'ATP synthase':1}
    bacteria = {'lipid metabolism': 1, 'Cofactor biosynthesis':1,
                'Metabolic intermediate biosynthesis': 1, 'ATP synthase':1}
    
    abundances = {'GCA_000275865.1_ASM27586v1': (6.4076277, archaea_co2),
                  'GCA_000013405.1_ASM1340v1': (2.3, bacteria), 
                  'GCA_000235565.1_ASM23556v1': (13.0963961, archaea_acetate), 
                  'GCA_000007985.2_ASM798v2': (0.9, bacteria), 
                  'GCA_000022125.1_ASM2212v1': (4.5, bacteria), 
                  'GCA_000970205.1_ASM97020v1': (19.6681987, archaea_acetate), 
                  'GCA_000012885.1_ASM1288v1': (2.2, bacteria), 
                  'GCA_000762265.1_ASM76226v1': (2.7, archaea_co2), 
                  'GCA_000016165.1_ASM1616v1': (4, bacteria), 
                  'GCA_000014965.1_ASM1496v1': (0.6473058, bacteria), 
                  'GCA_000014725.1_ASM1472v1': (3.0004899, bacteria), 
                  'GCA_000013445.1_ASM1344v1': (4.1752974, archaea_co2)}

    #genomes = glob.glob('../PostThesis/piptest/*.fna')
    rnaseqsim = RNASeqSim()
    for name in abundances.keys():
        genome = '../PostThesis/piptest/' + name + '_genomic.fna'
        output = '../PostThesis/piptest/new' + name
        print('to output:', output)
        path = pathlib.Path(output)
        path.mkdir(parents=True, exist_ok=True)
        info = rnaseqsim.join_information(genome, output)
        for factor in [1,3,1/6]:
            rnaseqsim.define_abundance(abundances[name], info, output, factor = factor)
            rnaseqsim.use_grinder(output, abundances[name][0], factor = factor)
            
    #samtools view -bS test1848/Assembly/mt1.sam | samtools sort - test1848/Assembly/mt1sorted
    
    perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome=../PostThesis/piptest/GCA_000275865.1_ASM27586v1_genomic.fna -out=../PostThesis/piptest/newGCA_000275865.1_ASM27586v1/fgs -complete=1 -train=./complete
    python MOSCA/rnaseqsim.py
    bowtie2 -f -x test1848/Assembly/contigs -1 polyester/sample_01_1.fasta -2 polyester/sample_01_2.fasta -S test1848/Analysis/sam1.sam
    bowtie2 -f -x test1848/Assembly/contigs -1 polyester/sample_02_1.fasta -2 polyester/sample_02_2.fasta -S test1848/Analysis/sam2.sam
    bowtie2 -f -x test1848/Assembly/contigs -1 polyester/sample_03_1.fasta -2 polyester/sample_03_2.fasta -S test1848/Analysis/sam3.sam
    
    htseq-count test1848/Analysis/sam1.sam test1848/Annotation/gff.gff.gff > sample1.readcounts
    htseq-count test1848/Analysis/sam2.sam test1848/Annotation/gff.gff.gff > sample2.readcounts
    htseq-count test1848/Analysis/sam3.sam test1848/Annotation/gff.gff.gff > sample3.readcounts
    
    
    #/home/jsequeira/biogrinder/script/grinder -reference_file fgs.ffn -abundance_file abundance1.config -total_reads 1000000 -mate_orientation FR -random_seed 13 -fastq_output 1 -qual_levels 30 10 -output_dir 1 -read_dist 251 -insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8
    
    '''
