# -*- coding: utf-8 -*-
"""
MOSCA's script for quantifying reads

By João Sequeira

Dec 2022
"""

import pandas as pd
import argparse
import multiprocessing
from mosca_tools import perform_alignment


def get_arguments():
    parser = argparse.ArgumentParser(description="MOSCA quantification")

    parser.add_argument("-o", "--output", help="Output directory (and input!).")
    parser.add_argument("-e", "--experiments", help="Filename of exps.")
    parser.add_argument(
        "-t", "--threads", default=multiprocessing.cpu_count() - 2, help="Number of threads to use [max available - 2]")
    args = parser.parse_args()
    args.output = args.output.rstrip('/')
    return args


def run():
    args = get_arguments()

    exps = pd.read_csv(args.experiments, sep='\t')

    for i in exps.index:
        if exps.iloc[i]['Data type'] == 'mrna':
            reference = f"{args.output}/Annotation/{exps.iloc[i]['Sample']}/fgs.ffn"
        elif exps.iloc[i]['Data type'] == 'dna':
            reference = f"{args.output}/Assembly/{exps.iloc[i]['Sample']}/contigs.fasta"
        else:
            continue
        perform_alignment(
            reference,
            [f"{args.output}/Preprocess/Trimmomatic/quality_trimmed_{exps.iloc[i]['Name']}_{fr}_paired.fq"
             for fr in ['forward', 'reverse']],
            f"{args.output}/Quantification/{exps.iloc[i]['Name']}", threads=args.threads)


if __name__ == '__main__':
    run()
