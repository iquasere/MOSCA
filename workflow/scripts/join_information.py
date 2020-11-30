import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import multiprocessing
from mosca_tools import timed_message, parse_blast, normalize_mg_readcounts_by_size, add_abundance, multi_sheet_excel, \
    normalize_readcounts, run_command

class Joiner:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's Protein and Entry reports")
        parser.add_argument("-e", "--experiments", type=str, required=True,
                            help="Experiments file")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads for reCOGnizer to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-if", "--input-format", type=str, default='tsv', choices=['tsv', 'excel'])
        parser.add_argument("-o", "--output", type=str, help="Output directory")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    def run(self):
        args = self.get_arguments()

        experiments = (pd.read_csv(args.experiments, sep='\t') if args.input_format == 'tsv' else
                       pd.read_excel(args.experiments))

        uniprotinfo = pd.read_csv('{}/Annotation/uniprotinfo.tsv'.format(args.output), sep='\t')
        taxonomy_columns = [col for col in uniprotinfo.columns if 'Taxonomic lineage' in col]
        functional_columns = ['COG general functional category', 'COG functional category',
                              'COG protein description', 'cog']

        mg_experiments = experiments[experiments["Data type"] == 'dna']

        sample2mgname = dict()
        for row in mg_experiments.iterrows():
            if row[1].loc['Sample'] in sample2mgname.keys():
                sample2mgname[row[1].loc['Sample']].append(row[1].loc['Name'])
            else:
                sample2mgname[row[1].loc['Sample']] = [row[1].loc['Name']]

        for sample in sample2mgname.keys():
            timed_message('Joining data for sample: {}'.format(sample))

            # Join BLAST and reCOGnizer outputs
            data = pd.merge(parse_blast('{}/Annotation/{}/aligned.blast'.format(args.output, sample)),
                            pd.read_excel('{}/Annotation/{}/protein2cog.xlsx'.format(args.output, sample)),
                            on='qseqid', how='left')
            data['sseqid'] = [ide.split('|')[1] if ide != '*' else ide for ide in data['sseqid']]
            data.columns = [column.replace('_x', ' (DIAMOND)').replace('_y', ' (reCOGnizer)') for column in
                            data.columns]  # after merging, BLAST columns will be repeated, and this makes explicit their origin
            data.columns = ['EC number (reCOGnizer)' if column == 'EC number' else column for column in data.columns]
            data = pd.merge(data, uniprotinfo, left_on='sseqid', right_on='Entry', how='left')
            data['Contig'] = [qseqid.split('_')[1] for qseqid in data['qseqid']]

            abundance_analysed = experiments[(experiments["Data type"] == 'dna') &
                                             (experiments["Sample"] == sample)]['Name'].tolist()
            expression_analysed = experiments[(experiments["Data type"] == 'mrna') &
                                                (experiments["Sample"] == sample)]['Name'].tolist()

            # MG quantification for each MG name of each Sample
            for mg_name in sample2mgname[sample]:
                # Normalization by contig size
                normalize_mg_readcounts_by_size(
                    '{}/Annotation/{}.readcounts'.format(args.output, mg_name),
                    '{}/Assembly/{}/contigs.fasta'.format(args.output, sample))

                data = add_abundance(data, '{}/Annotation/{}_normalized.readcounts'.format(
                                            args.output, mg_name), mg_name, origin_of_data='metagenomics',
                                     readcounts_has_tail=False)          # readcounts tail is removed in the normalization function

                for mt_name in expression_analysed:
                    data = add_abundance(data, '{}/Metatranscriptomics/{}.readcounts'.format(args.output, mt_name),
                                                mt_name, origin_of_data='metatranscriptomics')

            multi_sheet_excel('{}/MOSCA_Protein_Report.xlsx'.format(args.output), data, sheet_name=sample)

            print('Finding consensus COG for each Entry of Sample: {}'.format(sample))
            tqdm.pandas()
            cogs_df = data.groupby('Entry')['cog'].progress_apply(
                lambda x: x.value_counts().index[0] if len(x.value_counts().index) > 0
                else np.nan)
            cogs_df = cogs_df.reset_index()
            cogs_categories = data[functional_columns].drop_duplicates()

            # Aggregate information for each Entry, keep UniProt information, sum MG and MT or MP quantification
            data = data.groupby('Entry')[abundance_analysed + expression_analysed].sum().reset_index()
            data = pd.merge(data, uniprotinfo, on='Entry', how='left')
            data = pd.merge(data, cogs_df, on='Entry', how='left')
            data = pd.merge(data, cogs_categories, on='cog', how='left')
            data = data[uniprotinfo.columns.tolist() + functional_columns + abundance_analysed + expression_analysed]

            # MG normalization by sample and protein abundance
            data[abundance_analysed].to_csv(
                args.output + '/mg_preprocessed_readcounts.table', sep='\t', index=False)
            data = pd.concat([data, normalize_readcounts(args.output + '/mg_preprocessed_readcounts.table',
                    abundance_analysed, args.output + '/mg_preprocessed_normalization_factors.txt')[[
                    col + '_normalized' for col in abundance_analysed]]], axis=1)

            # MT normalization by sample and protein expression - normalization is repeated here because it's not stored from DESeq2 analysis
            data[expression_analysed].to_csv(args.output + '/expression_analysed_readcounts.table',
                                             sep='\t', index=False)
            data = pd.concat([data, normalize_readcounts(
                    args.output + '/expression_analysed_readcounts.table', expression_analysed,
                    args.output + '/expression_analysed_normalization_factors.txt')[[
                    col + '_normalized' for col in expression_analysed]]], axis = 1)

            # For each sample, write an Entry Report
            multi_sheet_excel('{}/MOSCA_Entry_Report.xlsx'.format(args.output),
                                     data, sheet_name=sample)

            data = pd.read_excel('{}/MOSCA_Entry_Report.xlsx'.format(args.output), sheet_name = sample + ' (1)')
            data[['Entry'] + expression_analysed].to_csv('{}/Metatranscriptomics/expression_matrix.tsv'.format(
                args.output), sep='\t', index=False)

            for mg_name in sample2mgname[sample]:
                # Draw the taxonomy krona plot
                data.groupby(taxonomy_columns)[mg_name].sum().reset_index()[[mg_name] + taxonomy_columns].to_csv(
                    '{}/{}_tax.tsv'.format(args.output, mg_name), sep='\t', index=False, header=False)
                run_command('ktImportText {0}/{1}_tax.tsv -o {0}/{1}_tax.html'.format(args.output, mg_name))

                # Draw the functional krona plot
                data.groupby(functional_columns)[mg_name].sum().reset_index()[[mg_name] + functional_columns].to_csv(
                    '{}/{}_fun.tsv'.format(args.output, mg_name), sep='\t', index=False, header=False)
                run_command('ktImportText {0}/{1}_fun.tsv -o {0}/{1}_fun.html'.format(args.output, mg_name))

if __name__ == '__main__':
    Joiner().run()