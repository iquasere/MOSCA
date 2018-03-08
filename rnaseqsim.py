# -*- coding: utf-8 -*-
"""
RNA-Seq simulator for differential expression

Created on Sun Jan 14 19:14:01 2018
"""

import analysing
import annotating
from diamond import DIAMOND
import subprocess
import pandas as pd

class RNASeqSim:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        
    def run(self, bashCommand):
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def use_fgs(self, genome, output):
        command = ('perl ../../../home/jsequeira/FGS/run_FragGeneScan.pl -genome='
                   + genome + ' -out=' + output + ' -complete=1 -train=./complete')
        self.run(command)
        
    def use_diamond(self, orfs, output):
        diamond = DIAMOND(threads = '6',
                          out = output,
                          query = orfs,
                          max_target_seqs = '1',
                          db = 'Databases/DIAMOND/UniProt/uniprot.dmnd')
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
        #self.use_fgs(genome, output + '/fgs')
        #self.use_diamond(output + '/fgs.faa', output + '/aligned.blast')
        blast = DIAMOND(out = output + '/aligned.blast').parse_result()
        blast = blast.drop_duplicates(subset = 'sseqid')
        ids = [ide.split('|')[1] if ide != '*' else ide for ide in blast['sseqid']]
        pathinfo = self.chunky(ids)
        blast.index = [ide.split('|')[1] if ide != '*' else ide for ide in blast['sseqid']]
        result = pd.concat([blast, pathinfo], axis=1, join='inner')
        result.index = result.qseqid
        result.to_csv(output + '/pathwayinfo', sep = '\t')
        return result[['qseqid', 'pathway','name']]

    #abundance is tuple (organism_presence, {pathway:expression})
    def define_abundance(self, abundance, info, output, factor = 1, base = 0):
        print(info.head())
        info['abundance'] = [0.01 for i in range(len(info))]
        for qseqid in info.index:
            for path, expression in abundance[1].items():
                if path in info.loc[qseqid]['pathway'] or 'ATP synthase' in info.loc[qseqid]['name']:
                    info.at[qseqid,'abundance'] = float(base + factor * expression * abundance[0])
        result = info[['abundance']]
        result.to_csv(output + '/abundance' + str(factor) + '.config', sep = '\t', header = False)
        
    def use_grinder(self, output, relative_abundance, factor = 1):
        command = ('/home/jsequeira/biogrinder/script/grinder -reference_file ' + output + '/fgs.ffn -abundance_file ' + 
        output + '/abundance' + str(factor) + '.config -total_reads ' + str(round(relative_abundance*100000)) + ' -mate_orientation FR -random_seed 13 -fastq_output 1 ' + 
        '-qual_levels 30 10 -output_dir ' + output + '/' + str(factor) + ' -read_dist 251 -insert_dist 2500 -mutation_dist ' + 
        'poly4 3e-3 3.3e-8')
        self.run(command)
