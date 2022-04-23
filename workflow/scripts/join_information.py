import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import multiprocessing
from mosca_tools import timed_message, parse_blast, normalize_mg_readcounts_by_size, add_abundance, multi_sheet_excel, \
    normalize_readcounts, run_command, parse_fasta


class Joiner:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's Protein and Entry reports")
        parser.add_argument("-e", "--experiments", type=str, required=True,
                            help="Experiments file")
        parser.add_argument("-t", "--threads", type=str, default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads for reCOGnizer to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-if", "--input-format", type=str, default='tsv', choices=['tsv', 'excel'])
        parser.add_argument("-o", "--output", type=str, help="Output directory")
        parser.add_argument("-nm", "--normalization-method", type=str, choices=["TMM", "RLE"],
                            help="Method for normalizing readcounts")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    def run(self):
        args = self.get_arguments()

        uniprotinfo = pd.read_csv(f'{args.output}/Annotation/uniprotinfo.tsv', sep='\t')

        sample2mgname = dict()
        for row in mg_experiments.iterrows():
            if row[1].loc['Sample'] in sample2mgname.keys():
                sample2mgname[row[1].loc['Sample']].append(row[1].loc['Name'])
            else:
                sample2mgname[row[1].loc['Sample']] = [row[1].loc['Name']]

        for sample in sample2mgname.keys():
            timed_message(f'Joining data for sample: {sample}')

            # new reCOGnizer exports results in EXCEL always
            recognizer_filename = f'{args.output}/Annotation/{sample}/reCOGnizer_results.xlsx'
            sheet_names = pd.ExcelFile(recognizer_filename).sheet_names
            cog_sheets = [sheet for sheet in sheet_names if 'COG' in sheet]
            cog_df = pd.read_excel(recognizer_filename, sheet_name=cog_sheets[0])
            for sheet in cog_sheets[1:]:
                cog_df = pd.concat([cog_df, pd.read_excel(recognizer_filename, sheet_name=sheet)])
            cog_df['cog'] = cog_df['DB ID']
            del cog_df['DB ID']

            # Join BLAST and reCOGnizer outputs
            data = pd.merge(parse_blast(f'{args.output}/Annotation/{sample}/aligned.blast'), cog_df,
                            on='qseqid', how='left')
            data['sseqid'] = data['sseqid_x']
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
                    f'{args.output}/Annotation/{mg_name}.readcounts',
                    f'{args.output}/Assembly/{sample}/contigs.fasta')

                data = add_abundance(data, f'{args.output}/Annotation/{mg_name}_normalized.readcounts', mg_name,
                                     origin_of_data='metagenomics', readcounts_has_tail=False)

            for mt_name in expression_analysed:
                data = add_abundance(data, f'{args.output}/Quantification/{mt_name}.readcounts', mt_name,
                                     origin_of_data='metatranscriptomics')

            multi_sheet_excel(f'{args.output}/MOSCA_Protein_Report.xlsx', data, sheet_name=sample)



            data[['Entry'] + expression_analysed].groupby('Entry')[expression_analysed].sum().reset_index().to_csv(
                f'{args.output}/Quantification/expression_matrix.tsv', sep='\t', index=False)

            for mg_name in sample2mgname[sample]:
                # Draw the taxonomy krona plot
                data.groupby(taxonomy_columns)[mg_name].sum().reset_index()[[mg_name] + taxonomy_columns].to_csv(
                    f'{args.output}/{mg_name}_tax.tsv', sep='\t', index=False, header=False)
                run_command('ktImportText {0}/{1}_tax.tsv -o {0}/{1}_tax.html'.format(args.output, mg_name))

                # Draw the functional krona plot
                data.groupby(functional_columns)[mg_name].sum().reset_index()[[mg_name] + functional_columns].to_csv(
                    f'{args.output}/{mg_name}_fun.tsv', sep='\t', index=False, header=False)
                run_command('ktImportText {0}/{1}_fun.tsv -o {0}/{1}_fun.html'.format(args.output, mg_name))


if __name__ == '__main__':
    Joiner().run()