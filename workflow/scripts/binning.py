# -*- coding: utf-8 -*-
"""
MOSCA's Binning package for clustering of 
contigs into Operational Taxonomic Units

By João Sequeira

Nov 2018
"""

from mosca_tools import run_command
import argparse
import multiprocessing

class Binner:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA binning")

        parser.add_argument("-c", "--contigs", type=str, required=True,
                            help="Filename of contigs")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use.")
        parser.add_argument("-o", "--output", type=str, help="Output directory"),
        parser.add_argument("-mset", "--markerset", type=str, default='40',
                            help="Set of marker genes to use (40 - all taxa; 107 - only bacteria)",
                            choices=['40', '107'])
        parser.add_argument("-s", "--sample", type=str, default='Sample',
                            help="Name of sample analysed")
        parser.add_argument("-r", "--reads", type=str, help="Filenames of reads")

        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        return args

    '''
    Input:
        contigs: FASTA file with contigs
        output: basename of output
        threads: number of threads to use by Maxbin
        mg1: name of forward reads file used in the assembly
        mg2: name of reverse reads file used in the assembly
        abundance: name of abundance file (format is contig\tabundance)
        marketset: either '107' marker genes present in >95% of bacteria, or
        '40' marker gene sets that are universal among bacteria and archaea. 
        '40' may be better suited for environment dominated by archaea; 
        however it tends to split genomes into more bins.
    Output:
        bins named basename + .n.fasta
        abundance of each contig in each sample (mg1 and mg2) named basename + .abund1/2
        abundace of each bin for both samples named basename + .abundance
        log of workflow named basename + .log
        markergenes used to compose the bins named basename + .marker
        contigs not included in any bin named basename + .noclass
    '''
    def run_maxbin(self, contigs, output, threads=8, reads=None, reads2=None,
                   abundance=None, markerset='40'):
        run_command('run_MaxBin.pl -contig {} -out {} -thread {} -markerset {}{}{}'.format(
            contigs, output, threads, markerset,' -reads1 {}'.format(reads) if reads else '',
            ' -reads2 {}'.format(reads2) if reads2 else ''))

    '''
    Input:
        bins_folder: str - foldername where the bins are
        output_directory: str - foldername where to store output
    Output:
        checkm.tsv is the table with completeness and contamination for each bin
    '''
    def run_checkm(self, bins_folder, threads='12'):
        run_command('checkm lineage_wf -x fasta -r --ali --nt -t {0} --pplacer_threads {0} {1} {1} --tab_table --file {1}/checkm.tsv'.format(
                threads, bins_folder))

    def run(self):
        args = self.get_arguments()

        files = args.reads.split(',')
        reads = files[0]
        if len(files) == 2:
            reads2 = files[1]

        sample = args.output.split('/')[-1]
        self.run_maxbin(args.contigs, '{}/{}'.format(args.output, sample), threads=args.threads,
                        reads=reads, reads2=reads2, markerset=args.markerset)
        self.run_checkm(args.output, threads=args.threads)

if __name__ == '__main__':
    Binner().run()