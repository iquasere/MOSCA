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
                            help="""Can be filename of single-end reads, comma separated-list of filenames of paired-end 
                            reads or filename of contigs file (this one requires the "--assembled option".""")
        parser.add_argument("-t", "--threads", type=str, default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-a", "--assembled", action="store_true", default=False,
                            help="If input is assembled reads.")
        parser.add_argument("-o", "--output", type=str, help="Output directory"),
        parser.add_argument("-em", "--error-model", type=str, default='illumina_5',
                            help="Error model for FastQ reads input",
                            choices=['sanger_5', 'sanger_10', '454_10', '454_30', 'illumina_5', 'illumina_10'])
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
        run_command(f"""run_FragGeneScan.pl -thread={threads} -genome={
            f'{file} -out={output} -complete=1 -train=./complete' if assembled else
            f'{file} -out={output} -complete=0 -train=./{error_model}'}""")

    def download_uniprot(self, out_dir):
        if not os.path.isfile(f'{out_dir}/uniprot.fasta'):
            for db in ['trembl', 'sprot']:
                run_command(
                    f'wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_'
                    f'{db}.fasta.gz')
            run_command('zcat uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz', output=f'{out_dir}/uniprot.fasta')
            run_command('rm uniprot_trembl.fasta.gz uniprot_sprot.fasta.gz')
        else:
            print(f'UniProt database found at: {out_dir}/uniprot.fasta')

    def generate_diamond_database(self, fasta, dmnd):
        run_command(f'diamond makedb --in {fasta} -d {dmnd}')

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
                self.generate_diamond_database(database, database.replace('fasta', 'dmnd'))
            else:
                print('DMND database was found. Using it')
            database = database.split('.fa')[0]
        elif database[-5:] != '.dmnd':
            print('Database must either be a FASTA (.fasta) or a DMND (.dmnd) file')
        else:
            if not os.path.isfile(database.replace('fasta', 'dmnd')):
                print('DMND database not found. Generating a new one')
                self.generate_diamond_database(database, database.replace('fasta', 'dmnd'))
            database = database.split('.dmnd')[0]

        run_command(
            f"diamond blastp --query {query} --out {aligned} --un {unaligned} --db {database} --outfmt 6 --unal 1 "
            f"--threads {threads} --max-target-seqs {max_target_seqs} -b {b} -c {c}")

    def run(self):
        args = self.get_arguments()

        pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

        self.gene_calling(args.input, f'{args.output}/fgs', threads=args.threads, assembled=args.assembled,
                          error_model=args.error_model)

        if args.download_uniprot:
            self.download_uniprot('/'.join(args.database.split('/')[:-1]))
            args.database = f"{'/'.join(args.database.split('/')[:-1])}/uniprot.fasta"

        # TODO - annotation has to be refined to retrieve better than hypothetical proteins
        (b, c) = self.b_n_c(argsb=args.block_size, argsc=args.index_chunks)
        self.run_diamond(f'{args.output}/fgs.faa', f'{args.output}/aligned.blast',
                         f'{args.output}/unaligned.fasta', args.database, threads=args.threads,
                         max_target_seqs=args.max_target_seqs, b=b, c=c)

if __name__ == '__main__':
    Annotater().run()