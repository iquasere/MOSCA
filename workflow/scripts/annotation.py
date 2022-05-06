# -*- coding: utf-8 -*-
"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

import argparse
import os
from multiprocessing import cpu_count
from mosca_tools import run_command
import pathlib


class Annotater:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA annotation")
        parser.add_argument(
            "-i", "--input", required=True,
            help='Can be filename of single-end reads, comma separated-list of filenames of paired-end reads or '
                 'filename of contigs/scaffolds file (the latter require the "--assembled" option.')
        parser.add_argument("-t", "--threads", default=cpu_count() - 2, help="Number of threads to use [available - 2]")
        parser.add_argument("--evalue", default=1e-3, help="Minimum e-value to report matches.")
        parser.add_argument(
            "-a", "--assembled", action="store_true", default=False, help="If input is assembled reads.")
        parser.add_argument("-o", "--output", help="Output directory")
        parser.add_argument(
            "-em", "--error-model", default='illumina_5', help="Error model for FastQ reads input",
            choices=['sanger_5', 'sanger_10', '454_10', '454_30', 'illumina_5', 'illumina_10', 'complete'])
        parser.add_argument("-db", "--database", help="Database for annotation")
        parser.add_argument(
            "--taxids", default=None, help="TaxIDs (comma-separated list) to retrieve reference proteomes of and use "
                                           "as the database for annotation.")
        parser.add_argument(
            "-mts", "--max-target-seqs", default=1, help="Number of identifications for each protein")
        parser.add_argument(
            "-b", "--block-size", default=None, help="Number of annotations to output per sequence inputed")
        parser.add_argument(
            "-c", "--index-chunks", default=None, help="Number of annotations to output per sequence inputed")
        parser.add_argument(
            "-cols", "--uniprot-columns", default=None, help="Columns to retrieve information from with UPIMAPI")
        parser.add_argument(
            "-dbs", "--uniprot-databases", default=None, help="Databases to cross-reference with UPIMAPI")
        parser.add_argument(
            "-rd", "--resources-directory", default=os.path.expanduser('~/resources'),
            help="Output directory for storing databases and other resources [~/resources]")


        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    def gene_calling(self, file, output, threads='12', assembled=True, error_model='illumina_10'):
        run_command(f"""run_FragGeneScan.pl -thread={threads} -genome={
            f'{file} -out={output} -complete=1 -train=./complete' if assembled else
            f'{file} -out={output} -complete=0 -train=./{error_model}'}""")

    def run_upimapi(
            self, query, output, rd='resources_directory', database='uniprot', evalue=0.001, threads=12,
            max_target_seqs=50, b=None, c=None, taxids=None, cols=None, dbs=None):
        run_command(
            f"upimapi.py -i {query} -o {output} -rd {rd} -db {database} --threads {threads} -mts {max_target_seqs} "
            f"--evalue {evalue}{f' --taxids {taxids}' if taxids is not None else ''}"
            f"{f' -b {b}' if b is not None else ''}{f' -c {c}' if c is not None else ''}"
            f"{f' --cols {cols}' if cols is not None else ''}{f' --dbs {dbs}' if dbs is not None else ''}")

    def run(self):
        args = self.get_arguments()

        pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

        self.gene_calling(
            args.input, f'{args.output}/fgs', threads=args.threads, assembled=args.assembled,
            error_model=args.error_model)

        self.run_upimapi(
            f'{args.output}/fgs.faa', args.output, rd=args.resources_directory, database=args.database, evalue=args.evalue, threads=args.threads,
            max_target_seqs=args.max_target_seqs, b=args.block_size, c=args.index_chunks, taxids=args.taxids)


if __name__ == '__main__':
    Annotater().run()
