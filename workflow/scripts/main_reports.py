# -*- coding: utf-8 -*-
"""
MOSCA's script for producing Protein and Entry reports

By João Sequeira

Dec 2022
"""

from pathlib import Path
from tqdm import tqdm
import argparse
import pandas as pd
import numpy as np
from mosca_tools import run_command, timed_message, multi_sheet_excel, normalize_mg_by_size, blast_cols, \
    normalize_readcounts


def get_arguments():
    parser = argparse.ArgumentParser(description="MOSCA main reports")

    parser.add_argument("-o", "--output", help="Output directory (and input!).")
    parser.add_argument("-e", "--experiments", help="Filename of exps.")

    args = parser.parse_args()
    args.output = args.output.rstrip('/')
    return args


def make_protein_report(out, exps):
    for sample in set(exps['Sample']):
        timed_message(f'Joining data for sample: {sample}')
        with open(f'{out}/Annotation/{sample}/fgs.faa') as f:
            lines = f.readlines()
        headers = [line.strip()[1:] for line in lines if line.startswith(">")]
        report = pd.DataFrame(headers, columns=["qseqid"])
        report = pd.merge(report, pd.read_csv(f'{out}/Annotation/{sample}/reCOGnizer_results.tsv', sep='\t'),
                          on='qseqid', how='left')
        report = report.groupby('qseqid')[report.columns.tolist()[1:]].first().reset_index()
        report = report[report['DB ID'].str.startswith('COG') == True].rename(columns={'DB ID': 'COG ID'})
        report = pd.merge(
            pd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t'), report, on='qseqid',
            how='outer')
        rename_cols = blast_cols + ['EC number']
        report = report.rename(columns={**{f'{col}_x': f'{col} (UPIMAPI)' for col in rename_cols},
                                        **{f'{col}_y': f'{col} (reCOGnizer)' for col in rename_cols}})
        report['Contig'] = report['qseqid'].apply(lambda x: x.split('_')[1])
        mg_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'dna')]['Name'].tolist()
        mt_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'mrna')]['Name'].tolist()
        mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()
        for mg_name in mg_names:
            readcounts = pd.read_csv(
                f'{out}/Quantification/{mg_name}.readcounts', sep='\t', header=None,
                names=['Contig', mg_name])
            normalize_mg_by_size(
                f'{out}/Quantification/{mg_name}.readcounts', f'{out}/Assembly/{sample}/contigs.fasta')
            norm_by_size = pd.read_csv(
                f'{out}/Quantification/{mg_name}_normalized.readcounts', sep='\t', header=None,
                names=['Contig', f'{mg_name} (Normalized by contig size)'])
            for counts in [readcounts, norm_by_size]:
                counts['Contig'] = counts['Contig'].apply(lambda x: x.split('_')[1])
            report = pd.merge(report, readcounts, on='Contig', how='left')
            report = pd.merge(report, norm_by_size, on='Contig', how='left')
        for mt_name in mt_names:
            readcounts = pd.read_csv(f'{out}/Quantification/{mt_name}.readcounts', sep='\t', header=None,
                                     names=['qseqid', mt_name])
            report = pd.merge(report, readcounts, on='qseqid', how='left')
        if len(mp_names) > 0:
            spectracounts = pd.read_csv(f'{out}/Metaproteomics/{sample}/spectracounts.tsv', sep='\t', header=None)
            report = pd.merge(report, spectracounts, on='qseqid', how='left')
        report[mg_names + mt_names + mp_names if len(mp_names) > 0 else []] = report[
            mg_names + mt_names + mp_names if len(mp_names) > 0 else []].fillna(value=0).astype(int)
        report[[f'{name} (Normalized by contig size)' for name in mg_names]].fillna(value=0, inplace=True)
        multi_sheet_excel(f'{out}/MOSCA_Protein_Report.xlsx', report, sheet_name=sample)
        report.to_csv(f'{out}/MOSCA_{sample}_Protein_Report.tsv', sep='\t', index=False)


def make_entry_report(protein_report, out, exps):
    for sample in set(exps['Sample']):
        timed_message(f'Organizing Entry level information for sample: {sample}')
        timed_message('Reading Protein Report')
        report = pd.read_excel(protein_report, sheet_name=sample)
        timed_message('Reading UPIMAPI report')
        upimapi_res = pd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t')
        uniprot_cols = [col for col in upimapi_res.columns if col not in blast_cols]
        taxonomy_columns = [col for col in upimapi_res.columns if 'Taxonomic lineage (' in col]
        functional_columns = ['General functional category', 'Functional category', 'Protein description', 'COG ID']
        if report['COG ID'].notnull().sum() > 0:
            tqdm.pandas(desc=timed_message(f'Finding consensus COG for each entry of sample: {sample}'))
            cogs_df = report.groupby('Entry')['COG ID'].progress_apply(
                lambda x: x.value_counts().index[0] if len(x.value_counts().index) > 0 else np.nan).reset_index()
            cogs_categories = report[functional_columns].drop_duplicates()
        else:
            timed_message('No COG information available')
            cogs_df = pd.DataFrame(columns=['Entry', 'COG ID'])
        mg_names = exps[(exps["Data type"] == 'dna') & (exps["Sample"] == sample)]['Name'].tolist()
        mt_names = exps[(exps["Data type"] == 'mrna') & (exps["Sample"] == sample)]['Name'].tolist()
        # Aggregate information for each Entry, keep UniProt information, sum MG and MT or MP quantification
        if len(mg_names) > 0:
            report.rename(
                columns={f'{mg_name} (Normalized by contig size)': mg_name for mg_name in mg_names}, inplace=True)
        report = report.groupby('Entry')[mg_names + mt_names].sum().reset_index()
        report = pd.merge(report, upimapi_res, on='Entry', how='left')
        report = pd.merge(report, cogs_df, on='Entry', how='left')
        if report['COG ID'].notnull().sum() > 0:
            report = pd.merge(report, cogs_categories, on='COG ID', how='left')
        else:
            report = pd.concat([report, pd.DataFrame(
                columns=['General functional category', 'Functional category', 'Protein description'],
                index=range(len(report)))], axis=1)
        report = report[uniprot_cols + functional_columns + mg_names + mt_names]
        timed_message('MG normalization by sample and protein abundance')
        if len(mg_names) > 0:
            report[mg_names].to_csv(f'{out}/Quantification/mg_preprocessed_readcounts.tsv', sep='\t', index=False)
            report = pd.concat([report, normalize_readcounts(
                f'{out}/Quantification/mg_preprocessed_readcounts.tsv', mg_names)[
                [f'{col}_normalized' for col in mg_names]]], axis=1)
        timed_message('MT normalization by sample and protein expression')
        if len(mt_names) > 0:
            report[mt_names].to_csv(f'{out}/Quantification/expression_analysed_readcounts.tsv', sep='\t',
                                    index=False)
            report = pd.concat([report, normalize_readcounts(
                f'{out}/Quantification/expression_analysed_readcounts.tsv', mt_names)[
                [f'{col}_normalized' for col in mt_names]]], axis=1)
        timed_message('Writing Entry Report')
        report = report.drop_duplicates()
        multi_sheet_excel(f'{out}/MOSCA_Entry_Report.xlsx', report, sheet_name=sample)
        report.to_csv(f'{out}/MOSCA_{sample}_Entry_Report.tsv', sep='\t', index=False)
        timed_message('Writting expression matrix')
        if len(mt_names) > 0:
            Path(f'{out}/Quantification/{sample}').mkdir(parents=True, exist_ok=True)
            report[['Entry'] + mt_names].groupby('Entry')[mt_names].sum().reset_index().to_csv(
                f'{out}/Quantification/{sample}/expression_matrix.tsv', sep='\t', index=False)
        if len(mg_names) == 0:
            mg_names = mt_names
        timed_message('Generating krona plots')
        for mg_name in mg_names:
            # Draw the taxonomy krona plot
            report.groupby(taxonomy_columns)[mg_name].sum().reset_index()[[mg_name] + taxonomy_columns].to_csv(
                f'{out}/{mg_name}_tax.tsv', sep='\t', index=False, header=False)
            run_command('ktImportText {0}/{1}_tax.tsv -o {0}/{1}_tax.html'.format(out, mg_name))
            # Draw the functional krona plot
            report.groupby(functional_columns)[mg_name].sum().reset_index()[[mg_name] + functional_columns].to_csv(
                f'{out}/{mg_name}_fun.tsv', sep='\t', index=False, header=False)
            run_command('ktImportText {0}/{1}_fun.tsv -o {0}/{1}_fun.html'.format(out, mg_name))


def run():
    args = get_arguments()
    exps = pd.read_csv(args.experiments, sep='\t')
    make_protein_report(args.output, exps),
    make_entry_report(f"{args.output}/MOSCA_Protein_Report.xlsx", args.output, exps)


if __name__ == '__main__':
    run()
