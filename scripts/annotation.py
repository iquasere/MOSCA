"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

from diamond import DIAMOND
from mosca_tools import MoscaTools
from uniprot_mapping import UniprotMapping
from progressbar import ProgressBar
import pandas as pd
import numpy as np
import os

mtools = MoscaTools()
upmap = UniprotMapping()

class Annotater:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    '''
    Input:
        file: name of input file to perform gene calling on
        output: basename of output files
        assembled: True if input is contigs, False if it are reads
        error_model: quality model to consider when input are reads
    Output:
        FragGeneScan output files will be produced with basename 'output + /fgs'
        If input is FASTQ reads (if assembled == False) a FASTA version will
        be produced in the same folder with the same name as the FASTQ file
        (but with .fasta instead of .fastq)
    '''
    def gene_calling(self, file, output, assembled = True, error_model = 'illumina_10'):
        bashCommand = 'run_FragGeneScan.pl -genome='
        if assembled:
            bashCommand += file + ' -out=' + output + '/fgs -complete=1 -train=./complete'
        else:
            mtools.fastq2fasta(file, file.replace('fastq', 'fasta'))            # fraggenescan only accepts FASTA input
            bashCommand += (file.replace('fastq', 'fasta') + ' -out=' + output +
                            '/fgs -complete=0 -train=./' + error_model)
        bashCommand += ' -thread=' + self.threads
        mtools.run_command(bashCommand)
    
    def annotation(self, max_target_seqs = '50'):
        diamond = DIAMOND(threads = self.threads,
                          db = self.db,
                          out = self.out_dir + '/aligned.blast',
                          query = self.out_dir + '/fgs.faa',
                          un = self.out_dir + '/unaligned.fasta',
                          unal = '1',
                          max_target_seqs = max_target_seqs)
        
        if self.db[-6:] == '.fasta' or self.db[-4:] == '.faa':
            print('FASTA database was inputed')
            if not os.path.isfile(self.db.replace('fasta','dmnd')):
                print('DMND database not found. Generating a new one')
                diamond.set_database(self.db, self.db.replace('fasta','dmnd'))
            else:
                print('DMND database was found. Using it')
        elif self.db[-5:] != '.dmnd':
            print('Database must either be a FASTA (.fasta) or a DMND (.dmnd) file')
            exit()
        diamond.db = self.db.split('.dmnd')[0] if '.dmnd' in self.db else self.db.split('.fasta')[0]
        
        diamond.run()
    
    def run(self):
        self.gene_calling(self.file, self.out_dir, self.assembled)
        #self.annotation()                                                      # annotation has to be refined to retrieved better than hypothetical proteins
        self.annotation(max_target_seqs = '1')

    '''
    UniProt regularly updates its databases. As we are working with MG here, many
    identifications will sometimes pertain to ORFs or pseudogenes that have been 
    wrongly predicted to code for proteins. It may also happen that the original
    authors decided to withdraw their published sequences.
    Input:
        uniprotinfo: name of UniProt info file
    Output:
        Rows in UniProt info file that lack a species identification will be 
        updated to include that new information
    '''
    def info_from_no_species(self, uniprotinfo):
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        missing_uniprotinfo = uniprotinfo[uniprotinfo['Taxonomic lineage (SPECIES)'].isnull()]
        ids = list(set([ide for ide in missing_uniprotinfo['Entry']]))
        new_uniprotinfo = self.get_uniprot_information(ids)
        for entry in new_uniprotinfo['Entry']:
            missing_uniprotinfo[missing_uniprotinfo['Entry'] == entry] = new_uniprotinfo[new_uniprotinfo['Entry'] == entry]
        new_uniprotinfo.to_csv(uniprotinfo, sep = '\t', index = False)
    
    '''
    Input:
        data: str - result from MOSCA analysis
        seqs_file: str - file to store FASTA sequences
        blast_file: str - file to store extra annotation
    Output:
        returns data with new protein names
    '''
    def further_annotation(self, data, temp_folder = 'temp',
                           dmnd_db = 'MOSCA/Databases/annotation_databases/uniprot.dmnd',
                           threads = '12', out_format = '6', all_info = True):
        '''
        noid_entries = data[(data['Protein names']=='Uncharacterized protein') & 
                       (data['COG general functional category']=='POORLY CHARACTERIZED')]['Entry']

        self.recursive_uniprot_fasta(temp_folder + '/temp.fasta', 
                                                   entries = noid_entries)
        
        mtools.run_command(('diamond blastp -q {} --db {} -p {} -f {} -o {}').format(
                temp_folder + '/temp.fasta', dmnd_db, threads, out_format, 
                temp_folder + '/temp_blast.tsv'))
        '''
        blast = mtools.parse_blast(temp_folder + '/temp_blast.tsv')
        '''
        noid_entries = [ide.split('|')[1] for ide in blast['sseqid']]
        
        uniprotinfo = self.get_uniprot_information(noid_entries)
        uniprotinfo.to_csv(temp_folder + '/up.info',index=False)
        '''
        uniprotinfo = pd.read_csv(temp_folder + '/up.info')
        uniprotinfo.drop_duplicates(inplace = True)
        
        relation = pd.DataFrame()
        relation['qseqid'] = [ide.split('|')[1] for ide in blast['qseqid']]
        relation['sseqid'] = [ide.split('|')[1] for ide in blast['sseqid']]
        relation = pd.merge(relation, uniprotinfo, left_on = 'sseqid', 
                            right_on = 'Entry', how = 'inner')
        relation = relation[(relation['Protein names'] != 'Deleted.') &
                            (relation['Protein names'] != 'Conserved protein') &
                            (~relation['Protein names'].str.contains('uncharacterized', case = False))]
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        data['Alternative entry'] = np.nan; data['Alternative protein name'] = np.nan
        pbar = ProgressBar()
        print('Searching for new annotations in the blast matches.')
        print('Proteins queried for new annotation:', len(set(relation['qseqid'])))
        print('Total possible identifications:', str(len(relation)))
        for query in pbar(relation['qseqid']):
            partial = relation[relation['qseqid'] == query]
            if all_info:
                data.loc[data['Entry'] == query, ['Alternative entry']] = '; '.join(partial['sseqid'])
                data.loc[data['Entry'] == query, ['Alternative protein name']] = '; '.join(partial['Protein names'])
            else:
                if len(partial) > 0:
                    original_taxonomy = data.loc[data['Entry'] == query].iloc[0][tax_columns]
                    max_proximity = 0; entry = np.nan; protein_name = np.nan
                    for i in range(len(partial)):
                        tax_score = self.score_taxonomic_proximity(original_taxonomy,
                            partial.iloc[i][tax_columns])
                        if tax_score > max_proximity:
                            max_proximity = tax_score
                            entry = partial.iloc[i]['sseqid']
                            protein_name = partial.iloc[i]['Protein names']
                    data.loc[data['Entry'] == query, ['Alternative entry']] = entry
                    data.loc[data['Entry'] == query, ['Alternative protein name']] = protein_name
        return data

    def score_taxonomic_proximity(self, tax1, tax2):
        proximity = 0
        for level in range(len(tax1)):
            if tax1[level] != tax2[level]:
                return proximity
            else:
                proximity += 1
        return proximity