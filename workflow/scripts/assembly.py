
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
from .mosca_tools import run_command, perform_alignment


class Assembler:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        if hasattr(self, 'memory') and self.assembler == 'metaspades':          # Megahit reads max memory in byte, metaspades reads in Gb
            self.memory /= 1.25e8
        
    def get_arguments(self):    
        parser = argparse.ArgumentParser(description="MOSCA assembly")
        parser.add_argument("-r", "--reads", type = str, required = True,
                            help="Reads files for assembly")
        parser.add_argument("-o", "--output", type = str, help = "Output directory")
        parser.add_argument("-t", "--threads", type = str, 
                            default = str(multiprocessing.cpu_count() - 2),
                            help = "Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-a","--assembler",type=str,choices=["metaspades","megahit"],
                            help="Tool for assembling the reads", default = "metaspades")
        parser.add_argument("-m","--memory",type=str, default = None,
                    help="Maximum memory (byte) available for assembly tools.")    
        args = parser.parse_args()
        
        args.output = args.output.rstrip('/')
        args.reads = args.reads.split(',')
            
        return args
        
    def run_assembler(self, reads, out_dir, assembler, threads = '12', 
                       memory = None):
        run_command('{} -o {} {} -t {}{}'.format(
            'metaspades.py' if assembler == 'metaspades' else 'megahit', out_dir, 
            '-1 {} -2 {}'.format(reads[0], reads[1]) if len(reads) == 2 else 
            '-{} {}'.format('s' if assembler == 'metaspades' else 'r', reads[0]), 
            threads, '-m {}'.format(memory) if memory else ''))
        
    def percentage_of_reads(self, file):
        handler = open(file)
        lines = handler.readlines()
        return lines[-1].split('%')[0]
    
    def run_metaquast(self, contigs, out_dir, threads = '12'):
        run_command('metaquast.py --threads {} --output-dir {} {}'.format(
            threads, out_dir, contigs))
     
    def run(self):
        args = self.get_arguments()
        
        # Assembly
        if args.assembler == 'megahit':                                         # snakemake workflow creates output directory, and megahit don't like it
            if os.path.isdir(args.output):
                shutil.rmtree(args.output)
                
        self.run_assembler(args.reads, args.output, args.assembler, 
                           threads = args.threads, memory = args.memory)
        if args.assembler == 'megahit':                                         # all contigs files are outputed the metaspades way
            os.rename('{}/final.contigs.fa'.format(args.output), 
                      '{}/contigs.fasta'.format(args.output))

        # Quality control
        pathlib.Path('{}/quality_control'.format(args.output)).mkdir(parents=True, exist_ok=True)

        self.run_metaquast('{}/contigs.fasta'.format(args.output), args.output)
        perform_alignment('{}/contigs.fasta'.format(args.output), args.reads, 
                          '{}/quality_control/alignment'.format(args.output), 
                          threads = args.threads)
        percentage_of_reads = self.percentage_of_reads('{}/quality_control/alignment.log'.format(args.output))
        
        if os.path.isfile('{}/quality_control/combined_reference/report.tsv'.format(args.output)):  #if metaquast finds references to the contigs, it will output results to different folders
            os.rename('{}/quality_control/combined_reference/report.tsv'.format(args.output), 
                      '{}/quality_control/report.tsv'.format(args.output))
        with open(args.output + '/quality_control/report.tsv', 'a') as f:
            f.write('Reads aligned (%)\t{}\n'.format(percentage_of_reads))
            
if __name__ == '__main__':
    Assembler().run()