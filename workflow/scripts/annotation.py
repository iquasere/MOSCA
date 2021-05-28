# -*- coding: utf-8 -*-
"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

import argparse
import multiprocessing
import numpy as np
import os
import pandas as pd
from mosca_tools import run_command
from progressbar import ProgressBar
import psutil
import pathlib


class Annotater:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA annotation")
        parser.add_argument("-i", "--input", type=str, required=True,
                            help="""Can be filename of single-end reads,
                            comma separated-list of filenames of paired-end reads
                            or filename of contigs file (this one requires the 
                            "--assembled option".""")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-a", "--assembled", action="store_true", default=False,
                            help="If input is assembled reads.")
        parser.add_argument("-o", "--output", type=str, help="Output directory"),
        parser.add_argument("-em", "--error-model", type=str, default='illumina_5',
                            help="Error model for FastQ reads input",
                            choices=['sanger_5', 'sanger_10', '454_10', '454_30',
                                     'illumina_5', 'illumina_10'])
        parser.add_argument("-db", "--database", type=str, help="Database for annotation")
        parser.add_argument("-mts", "--max-target-seqs", type=str, default=1,
                            help="Number of identifications for each protein")
        parser.add_argument("--download-uniprot", action="store_true", default=False,
                            help="Download uniprot database if FASTA DB doesn't exist")
        parser.add_argument("-b", "--block-size", default=None,
                            help="Number of annotations to output per sequence inputed")
        parser.add_argument("-c", "--index-chunks", default=None,
                            help="Number of annotations to output per sequence inputed")

        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        return args

    '''
    Input:
        file: name of input file to perform gene calling on
        output: basename of output files
        assembled: True if input is contigs, False if it are reads
        error_model: quality model to consider when input are reads
    Output:
        FragGeneScan output files will be produced with basename 'output'
        If input is FASTQ reads (if assembled == False) a FASTA version will
        be produced in the same folder with the same name as the FASTQ file
        (but with .fasta instead of .fastq)
    '''

    def gene_calling(self, file, output, threads='12', assembled=True, error_model='illumina_10'):
        run_command('run_FragGeneScan.pl -thread={} -genome={}'.format(threads,
                    '{} -out={} -complete=1 -train=./complete'.format(file, output) if assembled else
                    '{} -out={} -complete=0 -train=./{}'.format(file, output, error_model)))

    def download_uniprot(self, out_dir):
        if not os.path.isfile('{}/uniprot.fasta'.format(out_dir)):
            for db in ['trembl', 'sprot']:
                run_command(
                    'wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_{}.fasta.gz'.format(
                        db))
            run_command('zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz', output='{}/uniprot.fasta'.format(out_dir))
            run_command('rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz')
        else:
            print('UniProt database found at: {}/uniprot.fasta'.format(out_dir))

    def generate_diamond_database(self, fasta, dmnd):
        run_command('diamond makedb --in {} -d {}'.format(fasta, dmnd))

    def b_n_c(self, argsb, argsc):
        if argsb is not None:
            b = argsb
        else:
            b = psutil.virtual_memory().available / (1024.0 ** 3) / 20      # b = memory in Gb / 20
        if argsc is not None:
            return b, argsc
        if b > 3:
            return b, 1
        if b > 2:
            return b, 2
        if b > 1:
            return b, 3
        return b, 4

    def run_diamond(self, query, aligned, unaligned, database, threads=12, max_target_seqs=50, b=None, c=None):
        if database[-6:] == '.fasta' or database[-4:] == '.faa':
            print('FASTA database was inputed')
            if not os.path.isfile(database.replace('fasta', 'dmnd')):
                print('DMND database not found. Generating a new one')
                self.generate_diamond_database(database,
                                               database.replace('fasta', 'dmnd'))
            else:
                print('DMND database was found. Using it')
            database = database.split('.fa')[0]
        elif database[-5:] != '.dmnd':
            print('Database must either be a FASTA (.fasta) or a DMND (.dmnd) file')
        else:
            if not os.path.isfile(database.replace('fasta', 'dmnd')):
                print('DMND database not found. Generating a new one')
                self.generate_diamond_database(database,
                                               database.replace('fasta', 'dmnd'))
            database = database.split('.dmnd')[0]

        run_command(
            "diamond blastp --query {} --out {} --un {} --db {} --outfmt 6 --unal 1 --threads {} --max-target-seqs {} "
            "-b {} -c {}".format(query, aligned, unaligned, database, threads, max_target_seqs, b, c))

    def run(self):
        args = self.get_arguments()

        pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

        self.gene_calling(args.input, '{}/fgs'.format(args.output), threads=args.threads,
                          assembled=args.assembled, error_model=args.error_model)

        if args.download_uniprot:
            self.download_uniprot('/'.join(args.database.split('/')[:-1]))
            args.database = '{}/uniprot.fasta'.format('/'.join(args.database.split('/')[:-1]))

        # TODO - annotation has to be refined to retrieve better than hypothetical proteins
        (b, c) = self.b_n_c(argsb=args.block_size, argsc=args.index_chunks)
        self.run_diamond('{}/fgs.faa'.format(args.output), '{}/aligned.blast'.format(args.output),
                         '{}/unaligned.fasta'.format(args.output), args.database, threads=args.threads,
                         max_target_seqs=args.max_target_seqs, b=b, c=c)

    def further_annotation(self, data, temp_folder='temp',
                           dmnd_db='MOSCA/Databases/annotation_databases/uniprot.dmnd',
                           threads='12', out_format='6', all_info=True):
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
        uniprotinfo.drop_duplicates(inplace=True)

        relation = pd.DataFrame()
        relation['qseqid'] = [ide.split('|')[1] for ide in blast['qseqid']]
        relation['sseqid'] = [ide.split('|')[1] for ide in blast['sseqid']]
        relation = pd.merge(relation, uniprotinfo, left_on='sseqid',
                            right_on='Entry', how='inner')
        relation = relation[(relation['Protein names'] != 'Deleted.') &
                            (relation['Protein names'] != 'Conserved protein') &
                            (~relation['Protein names'].str.contains('uncharacterized', case=False))]
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)', 'Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)', 'Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)', 'Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        data['Alternative entry'] = np.nan;
        data['Alternative protein name'] = np.nan
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
                    max_proximity = 0
                    entry = np.nan
                    protein_name = np.nan
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


if __name__ == '__main__':
    Annotater().run()