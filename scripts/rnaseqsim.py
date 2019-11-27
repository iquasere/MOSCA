# -*- coding: utf-8 -*-
"""
RNA-Seq simulator for differential expression

Created on Sun Jan 14 19:14:01 2018
"""

from annotation import Annotater
from diamond import DIAMOND
import pandas as pd
from mosca_tools import MoscaTools
from annotation import Annotater
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
    def define_abundance_mt(self, excel, profiles, uniprotinfo, output_dir, 
                         factor = 1, base = 0):
        simulated = pd.read_excel(excel)
        simulated = simulated[simulated['Transcriptome link'].notnull()].reset_index()
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t', index_col = 0)
        handler = open('{}/abundance{}.config'.format(output_dir, str(factor)), 'w')
        print('Defining expressions for ' + str(len(simulated)) + ' organisms.')
        for i in range(len(simulated)):
            pbar = ProgressBar()
            print('Dealing with ' + simulated.iloc[i]['Species'])
            rna_file = 'SimulatedMGMT/rna/' + simulated.iloc[i]['Transcriptome link'].split('/')[-1].split('.gz')[0]
            abundance = simulated.iloc[i]['Abundance']
            profile = profiles[simulated.iloc[i]['Profile']]
            for rna in pbar(mtools.parse_fasta(rna_file).keys()):
                found = False
                ide = rna.split()[0]
                if ide not in uniprotinfo.index:
                    continue
                info = uniprotinfo.loc[ide]
                for path in profile['Pathway'].keys():
                    if path in str(info['Pathway']): 
                        handler.write(rna + '\t' + str(base + factor * abundance) + '\n')
                        found = True
                if not found:
                    if 'ATP synthase' in str(info['Protein names']):
                        handler.write(rna + '\t' + str(base + factor * abundance) + '\n')
                    else:
                        handler.write(rna + '\t' + str(base + abundance) + '\n')
        
    def use_grinder(self, reference_file, abundance_file, output_dir, coverage_fold = None, total_reads = None, seed = None):
        command = '{} -fastq_output 1 -qual_levels 30 10 -read_dist 151 -insert_dist 2500 -mutation_dist poly4 3e-3 3.3e-8 -mate_orientation FR -output_dir {} -reference_file {} -abundance_file {}{}{}{}'.format(
                'grinder', output_dir, reference_file, abundance_file, 
                   ' -total_reads ' + str(total_reads) if total_reads is not None else '',
                   ' -random_seed ' + seed if seed is not None else '',
                   ' -coverage_fold ' + coverage_fold if coverage_fold is not None else '',
                   
                   output_dir)
        mtools.run_command(command)
        
if __name__ == '__main__':
    import glob, pathlib
    
    data = pd.read_excel('SimulatedMGMT/simulated_taxa.xlsx')
    
    output_dir = 'SimulatedMGMT'
    
    rnaseqsimer = RNASeqSim()
    
    mtools = MoscaTools()
    
    # date is 14th November 2019
    '''
    # Download genomes and transcriptomes
    for i in range(len(data)):
        print(data.iloc[i]['Species'])
        if data.iloc[i]['Abundance'] > 0:
            print('Downloading ' + data.iloc[i]['Genome link'])
            mtools.run_command('wget {} -P {}/dna'.format(data.iloc[i]['Genome link'],
                           output_dir))
            if data.iloc[i]['Expression'] > 0:
                print('Downloading ' + data.iloc[i]['Transcriptome link'])
                mtools.run_command('wget {} -P {}/rna'.format(data.iloc[i]['Transcriptome link'],
                               output_dir))
    
    # Extract genomes and transcriptomes
    mtools.run_command('gunzip -v ' + ' '.join(glob.glob(output_dir + '/*/*.gz')))
    
    # Join all genomes into one file
    mg_files = glob.glob('SimulatedMGMT/dna/*.fa')
    #mtools.run_command('cat ' + ' '.join(mg_files), file = 'SimulatedMGMT/dna/genomes.fasta')
    # Join all transcriptomes into one file
    mt_files = glob.glob('SimulatedMGMT/rna/*.fa')
    #mtools.run_command('cat ' + ' '.join(mt_files), file = 'SimulatedMGMT/rna/transcriptomes.fasta')
    
    # Must retrieve data from UniProt manually, as I still don't know how to map from other DB than UniProt 
    # Data IDs are from Ensembl Genome Protein, as KZL88357 showed me
    # Get the IDs from genomes and transcriptomes
    print('Extracting IDs from ' + str(len(mg_files + mt_files)) + ' files.')
    ids = list()
    for file in (mg_files + mt_files):
        ids += mtools.parse_fasta(file).keys()
    
    ids = open('SimulatedMGMT/ids_for_uniprot_mapping.txt').read().split('\n')
    # Get the IDs to files to submit for UniProt
    print('Writing {} IDs to {}/ids_for_uniprot_mapping.txt'.format(str(len(ids)), output_dir))
    step = 10000
    j = 0
    pbar = ProgressBar()
    for i in pbar(range(0, len(ids), step)):
        open('{}/ids_for_uniprot_mapping{}.txt'.format(output_dir, str(j)), 'w').write('\n'.join(ids[i:i + step]))
        j += 1
    open('{}/ids_for_uniprot_mapping{}.txt'.format(output_dir, str(j)), 'w').write('\n'.join(ids[i:len(ids)]))
    
    uniprotinfo = pd.DataFrame()
    for file in glob.glob('SimulatedMGMT/*.tab'):
        df = pd.read_csv(file, sep = '\t')
        df.columns = ['Entry_name'] + df.columns.tolist()[1:]
        uniprotinfo = pd.concat([uniprotinfo, df])
        
    uniprotinfo.to_csv('SimulatedMGMT/uniprot_info.tsv', sep='\t', index=False)
    
    
    
    uniprotinfo = pd.read_csv('SimulatedMGMT/uniprot_info.tsv', sep='\t')
    
    profiles = {'archaea_co2': {
                        'Pathway':{
                                'One-carbon metabolism; methanogenesis from CO(2)': 1, 
                                'Cofactor biosynthesis': 1}, 
                        'Protein names':{
                                'ATP synthase': 1}
                        },
                'archaea_acetate': {
                        'Pathway':{
                                'One-carbon metabolism; methanogenesis from acetate': 1, 
                                'Cofactor biosynthesis': 1}, 
                        'Protein names':{
                                'ATP synthase': 1}
                        },
                'bacteria':{
                        'Pathway':{
                                'lipid metabolism': 1, 
                                'Cofactor biosynthesis': 1,
                                'Metabolic intermediate biosynthesis': 1}, 
                        'Protein names':{
                                'ATP synthase': 1}
                        }
                    }
                    
    for factor in [3,0.17,1]:
        rnaseqsimer.define_abundance_mt('SimulatedMGMT/simulated_taxa.xlsx',
           profiles, 'SimulatedMGMT/uniprot_info.tsv', 'SimulatedMGMT/rna', factor = factor)
    
    # grinder from 
    # perl Makefile.PL
    # make
    # sudo make install
    # cpan Bio::Perl
    
    mtools.run_command("sed -i 's/ /-/g' SimulatedMGMT/rna/*")
    
    seed = 13
    for factor in [1,3,0.17]:
        for letter in ['a','b','c']:
            pathlib.Path('{}/rna/{}/{}'.format(output_dir, letter, factor)).mkdir(parents=True, exist_ok=True)
            rnaseqsimer.use_grinder('SimulatedMGMT/rna/transcriptomes.fasta',
                'SimulatedMGMT/rna/abundance' + str(factor) + '.config',
                'SimulatedMGMT/rna/' + letter + '/' + str(factor), seed = str(seed),
                coverage_fold = '25')
            seed += 1
    '''
    '''
    with open(output_dir + '/dna/abundance.config', 'w') as f:
        for i in range(len(data)):
            if data.iloc[i]['Abundance'] > 0:
                file = '{}/dna/{}'.format(output_dir, data.iloc[i]['Genome link'].split('/')[-1].split('.gz')[0])
                seq_names = mtools.parse_fasta(file).keys()
                for seq_name in seq_names:
                    f.write('{}\t{}\n'.format(seq_name, str(float(data.iloc[i]['Abundance']) / len(seq_names))))
    
    import subprocess
    subprocess.check_output("sed -i 's/ /-/g' SimulatedMGMT/dna/*", shell=True)
    
    rnaseqsimer.use_grinder('SimulatedMGMT/dna/genomes.fasta',
                'SimulatedMGMT/dna/abundance.config', 'SimulatedMGMT/dna', 
                seed = '12', coverage_fold = '100')
    
    mtools.run_command('cat ' + ' '.join(glob.glob('SimulatedMGMT/dna/part*/grinder-reads.fastq')), 
                       file = 'SimulatedMGMT/dna/mg.fastq')
    mtools.divide_fq('SimulatedMGMT/dna/mg.fastq', 'SimulatedMGMT/dna/pretty_commune_R1.fastq', 'SimulatedMGMT/dna/pretty_commune_R2.fastq')
    
    for factor in [1,3,0.17]:
        for letter in ['a','b','c']:
            mtools.divide_fq('SimulatedMGMT/rna/{}/{}/grinder-reads.fastq'.format(letter, str(factor)),
                             'SimulatedMGMT/rna/{0}/{1}/rnaseq_{0}{1}reads_R1.fastq'.format(letter, str(factor)),
                             'SimulatedMGMT/rna/{0}/{1}/rnaseq_{0}{1}reads_R2.fastq'.format(letter, str(factor)))
    '''
    input_files = ' '.join(['SimulatedMGMT/dna/pretty_commune_R1.fastq,SimulatedMGMT/dna/pretty_commune_R2.fastq:' + 
          'SimulatedMGMT/rna/{0}/{1}/rnaseq_{0}{1}reads_R1.fastq,SimulatedMGMT/rna/{0}/{1}/rnaseq_{0}{1}reads_R2.fastq'.format(
                  letter, factor) for letter in ['a','b','c'] for factor in ['0.17','1','3']])
    
    conditions = 'c1,c2,c3,c1,c2,c3,c1,c2,c3'
    
    mtools.run_command('docker run -v /mnt/HDDStorage/jsequeira/:/input_data ' + 
        '-v /mnt/HDDStorage/jsequeira/SimulatedMGMT/:/MOSCA_analysis ' + 
        '-v /HDDStorage/jsequeira/MOSCA/Databases/annotation_databases/:/MOSCA/Databases/annotation_databases ' + 
        'iquasere/mosca -f {} -c {} -o SimulatedMGMT -t 14 -assstrat all'.format(input_files, conditions))