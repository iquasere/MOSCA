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
import pandas as pd
import shutil
import pathlib


class Binner:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA binning")

        parser.add_argument("-c", "--contigs", required=True, help="Filename of contigs")
        parser.add_argument(
            "-t", "--threads", default=multiprocessing.cpu_count() - 2, help="Number of threads to use.")
        parser.add_argument("-o", "--output", help="Output directory")
        parser.add_argument(
            "-mset", "--markerset", default='40', choices=['40', '107'],
            help="Set of marker genes to use (40 - all taxa; 107 - only bacteria)")
        parser.add_argument("-s", "--sample", default='Sample', help="Name of sample analysed")
        parser.add_argument("-r", "--reads", help="Filenames of reads")
        parser.add_argument("-ib", "--iterative-binning", action='store_true', default=False, help="Output directory")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    def run_maxbin(self, contigs, output, threads=8, reads=None, reads2=None, markerset='40', prob_threshold=0.9):
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
        output_folder = '/'.join(output.split('/')[:-1])  # remove the basename
        pathlib.Path(output_folder).mkdir(parents=True, exist_ok=True)
        run_command(
            f'run_MaxBin.pl -contig {contigs} -out {output} -thread {threads} -markerset {markerset} '
            f'-prob_threshold {prob_threshold}{f" -reads1 {reads}" if reads else ""}'
            f'{f" -reads2 {reads2}" if reads2 else ""}')

    def run_checkm(self, bins_folder, threads=12):
        run_command(
            f'checkm lineage_wf -x fasta -r --ali --nt -t {threads} --pplacer_threads {threads} {bins_folder} '
            f'{bins_folder} --tab_table --file {bins_folder}/checkm.tsv')

    def better_bin(self, table1, table2):
        table1 = pd.read_csv(table1, sep='\t')
        table2 = pd.read_csv(table2, sep='\t')

        hq_bins1 = ((table1.Completeness >= 90) & (table1.Contamination < 5)).sum()
        mq_bins1 = ((table1.Completeness < 90) & (table1.Completeness >= 50) & (table1.Contamination < 10)).sum()
        lq_bins1 = ((table1.Completeness < 50) & (table1.Contamination < 10)).sum()

        hq_bins2 = ((table2.Completeness >= 90) & (table2.Contamination < 5)).sum()
        mq_bins2 = ((table2.Completeness < 90) & (table2.Completeness >= 50) & (table2.Contamination < 10)).sum()
        lq_bins2 = ((table2.Completeness < 50) & (table2.Contamination < 10)).sum()

        if hq_bins1 > hq_bins2:
            return True
        if hq_bins1 < hq_bins2:
            return False
        if mq_bins1 > mq_bins2:
            return True
        if mq_bins1 < mq_bins2:
            return False
        if lq_bins1 > lq_bins2:
            return True
        if lq_bins1 < lq_bins2:
            return False
        return True

    def iterative_binning(self, contigs, output, threads=8, reads=None, reads2=None, markerset='40'):
        best_bin = 10
        sample = output.split('/')[-1]

        for prob_threshold in range(10, 100, 10):
            print(f'Probability threshold: {prob_threshold}')
            self.run_maxbin(contigs, f'{output}/{prob_threshold}/{sample}_{prob_threshold}', threads=threads,
                            reads=reads, reads2=reads2, markerset=markerset, prob_threshold=prob_threshold / 100)
            self.run_checkm(f'{output}/{prob_threshold}', threads=threads)

            if prob_threshold > 10:
                if self.better_bin(f'{output}/{prob_threshold}/checkm.tsv', f'{output}/{best_bin}/checkm.tsv'):
                    shutil.rmtree(f'{output}/{best_bin}')
                    print(f'Removed files for probability threshold: {best_bin} %')
                    best_bin = prob_threshold
                    print(f'New best probability threshold: {best_bin} %')
                else:
                    shutil.rmtree(f'{output}/{prob_threshold}')
                    print(f'Removed files for probability threshold: {prob_threshold} %')

        shutil.copyfile(f'{output}/{best_bin}/checkm.tsv', f'{output}/checkm.tsv')
        print(f'Best probability threshold: {best_bin} %')
        with open(f'{output}/result.txt', 'w') as f:
            f.write(f'Best probability threshold: {best_bin}')

    def run(self):
        args = self.get_arguments()

        files = args.reads.split(',')
        reads = files[0]
        if len(files) == 2:
            reads2 = files[1]

        sample = args.output.split('/')[-1]

        if args.iterative_binning:
            self.iterative_binning(
                args.contigs, args.output, threads=args.threads, reads=reads, reads2=reads2, markerset=args.markerset)
        else:
            self.run_maxbin(
                args.contigs, f'{args.output}/{sample}', threads=args.threads, reads=reads, reads2=reads2,
                markerset=args.markerset)
            self.run_checkm(args.output, threads=args.threads)


if __name__ == '__main__':
    Binner().run()
