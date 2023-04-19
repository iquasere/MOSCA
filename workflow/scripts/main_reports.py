"""
MOSCA's script for producing Protein and Entry reports

By João Sequeira

Dec 2022
"""

from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
from mosca_tools import run_command, timed_message, multi_sheet_excel, blast_cols


def make_protein_report(out, exps):
    for sample in set(exps['Sample']):
        timed_message(f'Joining data for sample: {sample}')
        with open(f'{out}/Annotation/{sample}/fgs.faa') as f:
            lines = f.readlines()
        headers = [line.strip()[1:] for line in lines if line.startswith(">")]
        report = pd.DataFrame(headers, columns=["qseqid"])
        cog_report = pd.read_csv(f'{out}/Annotation/{sample}/COG_report.tsv', sep='\t')
        cog_report = cog_report.groupby('qseqid')[cog_report.columns.tolist()[1:]].first().reset_index()
        report = pd.merge(report, cog_report, on='qseqid', how='left')
        report = report.groupby('qseqid')[report.columns.tolist()[1:]].first().reset_index()
        report = report[report['DB ID'].str.startswith('COG') == True].rename(columns={'DB ID': 'COG ID'})
        report = pd.merge(
            pd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t', low_memory=False),
            report, on='qseqid', how='outer')
        for col in report.columns:
            if '.' in col:      # some columns of UPIMAPI are repeated
                del report[col]
        rename_cols = blast_cols + ['EC number']
        report = report.rename(columns={**{f'{col}_x': f'{col} (UPIMAPI)' for col in rename_cols},
                                        **{f'{col}_y': f'{col} (reCOGnizer)' for col in rename_cols}})
        report['Contig'] = report['qseqid'].apply(lambda x: x.split('_')[1])
        mg_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'dna')]['Name'].tolist()
        mt_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'mrna')]['Name'].tolist()
        mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()
        for mg_name in mg_names:
            readcounts = pd.read_csv(
                f'{out}/Quantification/{mg_name}.readcounts', sep='\t', header=None, names=['Contig', mg_name])
            readcounts['Contig'] = readcounts['Contig'].apply(lambda x: x.split('_')[1])
            report = pd.merge(report, readcounts, on='Contig', how='left')
        for mt_name in mt_names:
            readcounts = pd.read_csv(
                f'{out}/Quantification/{mt_name}.readcounts', sep='\t', header=None, names=['qseqid', mt_name])
            report = pd.merge(report, readcounts, on='qseqid', how='left')
        if len(mt_names) > 0:
            report[['qseqid'] + mt_names].to_csv(f'{out}/Quantification/{sample}/mt.readcounts', sep='\t', index=False)
            print('Wrote MT readcounts.')
        if len(mp_names) > 0:
            spectracounts = pd.read_csv(f'{out}/Metaproteomics/{sample}/spectracounts.tsv', sep='\t')
            print(report.head())
            spectracounts.rename(columns={'Main Accession': 'qseqid'}, inplace=True)
            print(spectracounts.head())
            report = pd.merge(report, spectracounts, on='qseqid', how='left')
            report[['qseqid'] + mp_names][report[mp_names].isnull().sum(
                axis=1) < len(mp_names)].drop_duplicates().to_csv(
                f'{out}/Metaproteomics/{sample}/mp.spectracounts', sep='\t', index=False)
            print(report.head())
            print('Wrote MP spectracounts.')
        report[mg_names + mt_names + mp_names] = report[mg_names + mt_names + mp_names].fillna(value=0).astype(int)
        multi_sheet_excel(f'{out}/MOSCA_Protein_Report.xlsx', report, sheet_name=sample)
        report.to_csv(f'{out}/MOSCA_{sample}_Protein_Report.tsv', sep='\t', index=False)
        print(f'Wrote protein report for sample [{sample}].')


def make_entry_report(protein_report, out, exps):
    for sample in set(exps['Sample']):
        timed_message(f'Organizing Entry level information for sample: {sample}')
        timed_message('Reading Protein Report')
        report = pd.read_excel(protein_report, sheet_name=sample)
        timed_message('Reading UPIMAPI report')
        upimapi_res = pd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t')
        for col in upimapi_res.columns:
            if '.' in col:      # some columns of UPIMAPI are repeated
                del upimapi_res[col]
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
        mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()
        # Aggregate information for each Entry, keep UniProt information, sum MG and MT or MP quantification
        report = report.groupby('Entry')[mg_names + mt_names + mp_names].sum().reset_index()
        report = pd.merge(report, upimapi_res, on='Entry', how='left')
        report = pd.merge(report, cogs_df, on='Entry', how='left')
        if report['COG ID'].notnull().sum() > 0:
            report = pd.merge(report, cogs_categories, on='COG ID', how='left')
        else:
            report = pd.concat([report, pd.DataFrame(
                columns=['General functional category', 'Functional category', 'Protein description'],
                index=range(len(report)))], axis=1)
        report = report[uniprot_cols + functional_columns + mg_names + mt_names + mp_names]
        timed_message('Writing Entry Report')
        report = report.drop_duplicates()
        multi_sheet_excel(f'{out}/MOSCA_Entry_Report.xlsx', report, sheet_name=sample)
        report.to_csv(f'{out}/MOSCA_{sample}_Entry_Report.tsv', sep='\t', index=False)

        timed_message('Generating krona plots')
        tax_order = ['SUPERKINGDOM', 'KINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']
        print(taxonomy_columns)
        tax_cols = [
            f'Taxonomic lineage ({col})' for col in tax_order if f'Taxonomic lineage ({col})' in taxonomy_columns]
        print(tax_cols)
        for name in mg_names + mt_names + mp_names:
            Path(f'{out}/kronas').mkdir(parents=True, exist_ok=True)
            # Draw the taxonomy krona plot
            report.groupby(tax_cols)[name].sum().reset_index()[[name] + tax_cols].to_csv(
                f'{out}/kronas/{name}_tax.tsv', sep='\t', index=False, header=False)
            run_command('ktImportText {0}/kronas/{1}_tax.tsv -o {0}/kronas/{1}_tax.html'.format(out, name))
            # Draw the functional krona plot
            report.groupby(functional_columns)[name].sum().reset_index()[[name] + functional_columns].to_csv(
                f'{out}/kronas/{name}_fun.tsv', sep='\t', index=False, header=False)
            run_command('ktImportText {0}/kronas/{1}_fun.tsv -o {0}/kronas/{1}_fun.html'.format(out, name))


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')

    if snakemake.params.report == 'protein':
        make_protein_report(snakemake.params.output, exps)
    elif snakemake.params.report == 'entry':
        make_entry_report(f"{snakemake.params.output}/MOSCA_Protein_Report.xlsx", snakemake.params.output, exps)


if __name__ == '__main__':
    run()
