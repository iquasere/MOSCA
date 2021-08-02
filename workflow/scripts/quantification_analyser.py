# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package for retrieval of UniProt 
information and Differential Expression analysis

By João Sequeira

Sep 2017
"""

from mosca_tools import perform_alignment
import pandas as pd
import multiprocessing
import argparse
import os


class QuantificationAnalyser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA quantification analysis")
        parser.add_argument("-e", "--experiments", type=str, required=True,
                            help="Experiments file")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-if", "--input-format", type=str, default='tsv', choices=['tsv', 'excel'])
        parser.add_argument("-o", "--output", type=str, help="Output directory")

        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        return args

    '''
    input: 
        readcount_files: files from htseq-count with protein expression for each sample
        header: names of columns of final file / names of samples
        output: directory to output results
    output: 
        merged expression matrix name (output)
    '''
    def generate_expression_matrix(self, readcount_files, header, output):
        expression_matrix = pd.DataFrame()
        for file in readcount_files:
            df = pd.read_csv(file, sep='\t', index_col=0, header=None)
            df.columns = [file.split('/')[-1].rstrip('.readcounts')]
            expression_matrix = pd.merge(expression_matrix, df, how='outer',
                                         left_index=True, right_index=True)
        expression_matrix = expression_matrix[1:-5]  # remove non identified proteins (*) and the metrics at the end
        expression_matrix = expression_matrix.fillna(value=0).astype(
            int)  # set not identified proteins expression to 0, and convert floats, since DeSEQ2 only accepts integers
        expression_matrix.index.name = 'geneid'
        expression_matrix.columns = header
        expression_matrix.to_csv(output, sep='\t')

    def run(self):
        args = self.get_arguments()

        experiments = (pd.read_csv(args.experiments, sep='\t') if args.input_format == 'tsv' else
                       pd.read_excel(args.experiments))
        mt_experiments = experiments[experiments['Data type'] == 'mrna']
        for i in experiments.index:
            if experiments.iloc[i]['Data type'] == 'mrna':
                reference = f"{args.output}/Annotation/{experiments.iloc[i]['Sample']}/fgs.ffn"
            elif experiments.iloc[i]['Data type'] == 'dna':
                reference = f"{args.output}/Assembly/{experiments.iloc[i]['Sample']}/contigs.fasta"
            elif experiments.iloc[i]['Data type'] == 'protein':
                continue
            else:
                print('A data type MOSCA can yet not handle!')
                continue

            print(f"Data type is: {experiments.iloc[i]['Data type']}")

            if not os.path.isfile(f"{args.output}/Quantification/{experiments.iloc[i]['Name']}.readcounts"):
                print(f"{args.output}/Quantification/{experiments.iloc[i]['Name']}.readcounts not found! Generating it")
                perform_alignment(
                    reference,
                    [f"{args.output}/Preprocess/Trimmomatic/"
                     f"quality_trimmed_{experiments.iloc[i]['Name']}_{fr}_paired.fq" for fr in ['forward', 'reverse']],
                    f"{args.output}/Quantification/{experiments.iloc[i]['Name']}",
                    threads=args.threads)
            else:
                print(f"{args.output}/Quantification/{experiments.iloc[i]['Name']}.readcounts exists.")

        self.generate_expression_matrix(
            [f"{args.output}/Quantification/{mt_name}.readcounts" for mt_name in mt_experiments['Name']],
            mt_experiments['Name'].tolist(),
            f"{args.output}/Quantification/expression_matrix.tsv")


if __name__ == '__main__':
    QuantificationAnalyser().run()
