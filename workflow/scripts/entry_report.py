"""
MOSCA's script for producing Protein and Entry reports

By João Sequeira

Dec 2022
"""

from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import csv
from mosca_tools import run_command, timed_message, multi_sheet_excel

functional_columns = [
    'General functional category', 'Functional category', 'Protein description', 'COG ID', 'EC number (reCOGnizer)']


def estimate_cog_for_entries(e_report):
    if e_report['COG ID'].notnull().sum() == 0:
        print('No COG information available.')
        return pd.DataFrame()
    tqdm.pandas(desc=timed_message(f'Finding consensus COG for each entry.'), ascii=' >=')
    cogs_df = e_report.groupby('Entry')['COG ID'].progress_apply(
        lambda x: x.value_counts().index[0] if len(x.value_counts().index) > 0 else np.nan).reset_index()
    cog_categories = e_report[functional_columns].drop_duplicates()
    duplicated_cogs = cog_categories[cog_categories['COG ID'].duplicated()]['COG ID'].unique()
    categories_df = pd.concat([
        cog_categories[~cog_categories['COG ID'].isin(duplicated_cogs)],
        cog_categories[cog_categories['COG ID'].isin(duplicated_cogs)][~(
                cog_categories['General functional category'].isnull() & cog_categories['COG ID'].notnull())]])
    fcols = functional_columns.copy()
    fcols.remove('COG ID')
    categories_df = categories_df.groupby('COG ID')[fcols].first().reset_index()
    cogs_df = pd.merge(cogs_df, categories_df, on='COG ID', how='left')
    cogs_df = cogs_df.groupby('Entry')[functional_columns].first().reset_index()    # Keep only one COG per Entry
    return cogs_df


def join_normalized_matrices(mg_names, mt_names, mp_names, out):
    counts = pd.DataFrame()
    if len(mg_names) > 0:
        counts = pd.merge(counts, pd.read_csv(
            f'{out}/Quantification/mg_normalized.tsv', sep='\t'), left_index=True, right_index=True, how='outer')
    if len(mt_names) > 0:
        counts = pd.merge(counts, pd.read_csv(
            f'{out}/Quantification/mt_normalized.tsv', sep='\t'), left_index=True, right_index=True, how='outer')
    if len(mp_names) > 0:
        counts = pd.merge(counts, pd.read_csv(
            f'{out}/Metaproteomics/mp_normalized.tsv', sep='\t'), left_index=True, right_index=True, how='outer')
    counts.rename(columns={col: f'{col}_normalized' for col in counts.columns.tolist()}, inplace=True)
    if 'Entry' not in counts.columns.tolist():
        counts = counts.reset_index().rename(columns={'index': 'Entry'})
    return counts


def write_kronas(report, out, mg_names, mt_names, mp_names, tax_cols):
    timed_message('Generating krona plots')
    #tax_order = ['SUPERKINGDOM', 'KINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']
    #tax_cols = [f'Taxonomic lineage ({col})' for col in tax_order if f'Taxonomic lineage ({col})' in taxonomy_columns]
    rep = report.copy()
    rep[tax_cols + functional_columns] = rep[tax_cols + functional_columns].fillna(value='Unknown')
    for name in mg_names + mt_names + mp_names:
        Path(f'{out}/kronas').mkdir(parents=True, exist_ok=True)
        # Draw the taxonomy krona plot
        rep_name = rep.groupby(tax_cols)[name].sum().reset_index()[[name] + tax_cols]
        rep_name = rep_name[rep_name[name] > 0]
        rep_name.to_csv(f'{out}/kronas/{name}_tax.tsv', sep='\t', index=False, header=False)
        run_command(
            'ktImportText {0}/kronas/{1}_tax.tsv -o {0}/kronas/{1}_tax.html'.format(out, name),
            print_message=False)
        # Draw the functional krona plot
        rep_name = rep.groupby(functional_columns)[name].sum().reset_index()[[name] + functional_columns]
        rep_name = rep_name[rep_name[name] > 0]
        rep_name.to_csv(f'{out}/kronas/{name}_fun.tsv', sep='\t', index=False, header=False)
        run_command(
            'ktImportText {0}/kronas/{1}_fun.tsv -o {0}/kronas/{1}_fun.html'.format(out, name),
            print_message=False)


def get_lowest_otu(entry, tax_cols):
    i = 0
    for col in tax_cols[::-1]:
        if not pd.isna(entry[col]):
            if i == 0 or 'Candidatus' in entry[col] or 'Other' in entry[col]:
                return entry[col]
            return f'Other {entry[col]}'
        i += 1
    return np.nan


def make_entry_report(out, exps):
    mg_names = exps[(exps["Data type"] == 'dna')]['Name'].tolist()
    mt_names = exps[(exps["Data type"] == 'mrna')]['Name'].tolist()
    mp_names = exps[(exps['Data type'] == 'protein')]['Name'].tolist()
    with open(f"{out}/Annotation/{exps['Sample'].tolist()[0]}/uniprotinfo.tsv", 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        up_cols = next(reader)
    all_info = pd.DataFrame()
    for sample in set(exps['Sample']):
        timed_message(f'Reading Protein Report for sample: {sample}.')
        all_info = pd.concat([all_info, pd.read_csv(
            f'{out}/MOSCA_{sample}_Protein_Report.tsv', sep='\t', low_memory=False)])
    if 'EC number' in up_cols:
        up_cols[up_cols.index('EC number')] = 'EC number (UPIMAPI)'
    entry_report = all_info[up_cols].drop_duplicates()
    entry_report = pd.merge(entry_report, estimate_cog_for_entries(all_info), on='Entry', how='left')
    timed_message('Adding quantification at the entry level.')
    entry_report = pd.merge(
        entry_report, all_info.groupby('Entry')[mg_names + mt_names + mp_names].sum().reset_index(), on='Entry',
        how='left')
    timed_message('Adding normalized matrices.')
    entry_report = pd.merge(
        entry_report, join_normalized_matrices(mg_names, mt_names, mp_names, out), on='Entry', how='left')
    timed_message('Writing Entry Report.')
    entry_report.to_csv(f'{out}/MOSCA_Entry_Report.tsv', sep='\t', index=False)
    multi_sheet_excel(f'{out}/MOSCA_Entry_Report.xlsx', entry_report, sheet_name='Entry Report')
    taxonomy_columns = [col for col in up_cols if 'Taxonomic lineage (' in col and 'IDS' not in col.upper()]
    timed_message('Writing KEGGcharter input.')
    entry_report.to_csv(f'{out}/keggcharter_input.tsv', sep='\t', index=False)
    timed_message('Writing Krona plots.')
    for i in range(len(taxonomy_columns)):
        entry_report[taxonomy_columns[i]] = entry_report.apply(
            lambda x: get_lowest_otu(x, taxonomy_columns[:i + 1]), axis=1)
    entry_report['Taxonomic lineage (SPECIES)'] = entry_report.apply(
        lambda x: get_lowest_otu(x, taxonomy_columns), axis=1)
    write_kronas(
        entry_report, out, [f'{name}_normalized' for name in mg_names],
        [f'{name}_normalized' for name in mt_names],
        [f'{name}_normalized' for name in mp_names], taxonomy_columns)


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')
    make_entry_report(snakemake.params.output, exps)


if __name__ == '__main__':
    run()
