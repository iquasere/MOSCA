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

        experiments = pd.read_csv(args.experiments, sep='\t')
        mt_experiments = experiments[experiments['Data type'] == 'mrna']

        for i in mt_experiments.index:
            if experiments.iloc[i]['Data type'] == 'mrna':
                attribute = 'Name'
                folder = 'Metatranscriptomics'
            else:
                attribute = 'gene_id'
                folder = 'Annotation'

            if not os.path.isfile("{}/{}/{}.readcounts".format(args.output, folder, experiments.iloc[i]['Name'])):
                perform_alignment('{}/Assembly/{}/contigs.fasta'.format(args.output, experiments.iloc[i]['Sample']),
                              ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_{}_paired.fq'.format(
                                args.output, experiments.iloc[i]['Name'], fr) for fr in ['forward', 'reverse']],
                              '{}/{}/{}'.format(args.output, folder, experiments.iloc[i]['Name']),
                              blast=('{}/Annotation/{}/aligned.blast'.format(args.output, experiments.iloc[i]['Sample'])
                                     if attribute == 'Name' else None), threads=args.threads, attribute=attribute)

        self.generate_expression_matrix(["{}/Metatranscriptomics/{}.readcounts".format(args.output, mt_name) for mt_name in
                                         mt_experiments['Name']], mt_experiments['Name'].tolist(),
                                        '{}/Metatranscriptomics/expression_matrix.tsv'.format(args.output))

if __name__ == '__main__':
    QuantificationAnalyser().run()