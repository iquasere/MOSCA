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
from mosca_tools import run_command, run_pipe_command, perform_alignment
import psutil


class Assembler:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs


    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA assembly")
        parser.add_argument("-r", "--reads", type=str, required=True,
                            help="Reads files for assembly")
        parser.add_argument("-o", "--output", type=str, help="Output directory")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-a", "--assembler", type=str, choices=["metaspades", "megahit", "rnaspades"],
                            help="Tool for assembling the reads", default="metaspades")
        # default memory is a third of total available memory
        parser.add_argument("-m", "--memory", default=psutil.virtual_memory().available / (1024.0 ** 3) / 3, type=float,
                            help="Maximum memory (Gb) available for assembly tools.")
        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        args.reads = args.reads.split(',')

        if args.assembler == 'megahit':     # Megahit reads max memory in byte
            args.memory *= 10e9

        return args

    def run_assembler(self, reads, out_dir, assembler, threads='12', memory=None):
        run_command(f"{assembler if assembler == 'megahit' else f'{assembler}.py'} -o {out_dir} -t {threads} "
                    f"""{f'-1 {reads[0]} -2 {reads[1]}' if len(reads) == 2 else 
                        f'-{"r" if assembler == "megahit" else "s"} {reads[0]}'} """
                    f"-m {round(memory) if memory else ''}")

    def percentage_of_reads(self, file):
        handler = open(file)
        lines = handler.readlines()
        return lines[-1].split('%')[0]

    def run_metaquast(self, contigs, out_dir, threads='12'):
        run_command(f'metaquast.py --threads {threads} --output-dir {out_dir} --max-ref-number 0 {contigs}')

    def run(self):
        args = self.get_arguments()

        # Assembly
        if args.assembler == 'megahit':  # snakemake workflow creates output directory, and megahit don't like it
            if os.path.isdir(args.output):
                shutil.rmtree(args.output)

        self.run_assembler(args.reads, args.output, args.assembler, threads=args.threads, memory=args.memory)

        if args.assembler == 'megahit':  # all contigs files are outputed the metaspades way
            run_pipe_command(f"awk \'{{print $1}}\' {args.output}/final.contigs.fa",
                             # k141_714 flag=1 multi=1.0000 len=369 -> k141_714
                             output=f'{args.output}/contigs.fasta')

        # Quality control
        pathlib.Path(f'{args.output}/quality_control').mkdir(parents=True, exist_ok=True)

        self.run_metaquast(f'{args.output}/contigs.fasta', f'{args.output}/quality_control')
        perform_alignment(f'{args.output}/contigs.fasta', args.reads, f'{args.output}/quality_control/alignment',
                          threads=args.threads)
        percentage_of_reads = self.percentage_of_reads(f'{args.output}/quality_control/alignment.log')

        if os.path.isfile(f'{args.output}/quality_control/combined_reference/report.tsv'
                        ):  # if metaquast finds references to the contigs, it will output results to different folders
            shutil.copyfile(f'{args.output}/quality_control/combined_reference/report.tsv',
                            f'{args.output}/quality_control/report.tsv')
        with open(args.output + '/quality_control/report.tsv', 'a') as f:
            f.write(f'Reads aligned (%)\t{percentage_of_reads}\n')


if __name__ == '__main__':
    Assembler().run()
