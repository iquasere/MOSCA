"""
MOSCA's Assembly package for performing Assembly with MetaSPAdes
and Megahit and Quality Control analysis of the resulting contigs

By João Sequeira

Jun 2017
"""

import argparse
import multiprocessing
import os
import pathlib
import shutil
from glob import glob

from mosca_tools import run_command, run_pipe_command, perform_alignment
import psutil


class Assembler:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA assembly")
        parser.add_argument("-r", "--reads", required=True, help="Reads files for assembly")
        parser.add_argument("-o", "--output", help="Output directory")
        parser.add_argument(
            "-t", "--threads", default=multiprocessing.cpu_count() - 2,
            help="Number of threads to use [max available - 2]")
        parser.add_argument(
            "-a", "--assembler", choices=["metaspades", "megahit", "rnaspades", "trinity"],
            help="Tool for assembling the reads [metaspades]", default="metaspades")
        parser.add_argument(
            "-m", "--memory", default=psutil.virtual_memory().available / (1024.0 ** 3) / 3, type=float,
            help="Maximum memory (Gb) available for assembly tools [max available memory / 3]")
        parser.add_argument("-rl", "--read-length", default=150, help="Average length of reads [150]")
        parser.add_argument("-is", "--insert-size", default=2500, help="Average insert size [2500]")
        parser.add_argument(
            "-stdi", "--std-insert", default=10, help="Standard deviation of insert size [2500]")
        parser.add_argument(
            "-bq", "--base-qual", default='phred33', choices=['phred33', 'phred64'],
            help="Quality format of base call of reads files")
        parser.add_argument(
            "-mrn", "--max-ref-number", type=int, default=50, help="Maximum references to be downloaded by MetaQUAST")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        args.reads = args.reads.split(',')
        if args.assembler == 'megahit':     # Megahit reads max memory in byte
            args.memory *= 10e9
        return args

    def run_assembler_mg(self, reads, out_dir, assembler, threads='12', memory=None):
        run_command(
            f"{assembler if assembler == 'megahit' else f'{assembler}.py'} -o {out_dir} -t {threads} "
            f"""{f'-1 {reads[0]} -2 {reads[1]}' if len(reads) == 2 else 
            f'-{"r" if assembler == "megahit" else "s"} {reads[0]}'} """
            f"-m {round(memory) if memory else ''}")

    def run_assembler_mt(self, reads, out_dir, threads, memory=None):
        reads = f'--left {reads[0]} --right {reads[1]}' if len(reads) == 2 else f'--single {reads[0]}'
        run_command(
            f"Trinity --seqType fq {reads} --CPU {threads} --max_memory {int(memory)}G --output {out_dir}/trinity")

    def percentage_of_reads(self, file):
        handler = open(file)
        lines = handler.readlines()
        return lines[-1].split('%')[0]

    def run_metaquast(self, contigs, out_dir, threads='12', max_ref_number=0):
        run_command(
            f'metaquast.py --threads {threads} --output-dir {out_dir} --max-ref-number {max_ref_number} {contigs}')

    def close_gaps(
            self, scaffolds, contigs, output_basename, read1, read2, read_length=150, insert_size=1500, std_insert=10,
            threads=14, base_qual='phred33'):
        run_command(
            f'gmcloser --target_scaf {scaffolds} --query_seq {contigs} --prefix_out {output_basename} '
            f'-r {read1} {read2} --read_len {read_length} --insert {insert_size} --sd_insert {std_insert} '
            f'--thread {threads} --base_qual {base_qual}')

    def run(self):
        args = self.get_arguments()

        # Assembly
        if args.assembler == 'megahit':  # snakemake workflow creates output directory, and megahit don't like it
            if os.path.isdir(args.output):
                shutil.rmtree(args.output)

        if args.assembler != 'trinity':
            self.run_assembler_mg(args.reads, args.output, args.assembler, threads=args.threads, memory=args.memory)
        else:
            self.run_assembler_mt(args.reads, args.output, threads=args.threads, memory=args.memory)
        print(f'Assembler is: {args.assembler}')
        if args.assembler == 'megahit':  # all contigs files are outputed the metaspades way
            run_pipe_command(f"awk \'{{print $1}}\' {args.output}/final.contigs.fa",
                             # k141_714 flag=1 multi=1.0000 len=369 -> k141_714
                             output=f'{args.output}/contigs.fasta')
            shutil.copyfile(f'{args.output}/contigs.fasta', f'{args.output}/scaffolds.fasta'
                            )   # TODO - put SOAPdenovo producing scaffolds from Megahit
        elif args.assembler == 'trinity':
            run_pipe_command(f"awk \'{{print $1}}\' {args.output}/trinity/Trinity.fasta",
                             # >TRINITY_DN19419_c0_g1_i1 len=248 path=[0:0-247] -> TRINITY_DN19419_c0_g1_i1
                             output=f'{args.output}/contigs.fasta')
            shutil.copyfile(f'{args.output}/contigs.fasta', f'{args.output}/scaffolds.fasta')

        self.close_gaps(f'{args.output}/scaffolds.fasta', f'{args.output}/contigs.fasta', f'{args.output}/gap_close',
                        args.reads[0], args.reads[1], read_length=args.read_length, insert_size=args.insert_size,
                        std_insert=args.std_insert, threads=args.threads, base_qual=args.base_qual)

        # Quality control
        pathlib.Path(f'{args.output}/quality_control').mkdir(parents=True, exist_ok=True)

        self.run_metaquast(f'{args.output}/contigs.fasta', f'{args.output}/quality_control')
        perform_alignment(f'{args.output}/contigs.fasta', args.reads, f'{args.output}/quality_control/alignment',
                          threads=args.threads)
        percentage_of_reads = self.percentage_of_reads(f'{args.output}/quality_control/alignment.log')

        if os.path.isfile(f'{args.output}/quality_control/combined_reference/report.tsv'
                          ):  # if metaquast finds references to contigs, it will output results to different folders
            shutil.copyfile(f'{args.output}/quality_control/combined_reference/report.tsv',
                            f'{args.output}/quality_control/report.tsv')
        with open(args.output + '/quality_control/report.tsv', 'a') as f:
            f.write(f'Reads aligned (%)\t{percentage_of_reads}\n')


if __name__ == '__main__':
    Assembler().run()
