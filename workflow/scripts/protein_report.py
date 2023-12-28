"""
MOSCA's script for producing Protein report

By João Sequeira

Dec 2022
"""

import shutil
import pandas as pd
from mosca_tools import timed_message, blast_cols

functional_columns = [
    'General functional category', 'Functional category', 'Protein description', 'COG ID', 'EC number (reCOGnizer)']


def make_protein_report(out, exps, sample, mg_preport, mt_preport, mp_preport, de_input):
    timed_message(f'Joining data for sample: {sample}.')
    with open(f'{out}/Annotation/{sample}/fgs.faa') as f:
        lines = f.readlines()
    headers = [line.strip()[1:] for line in lines if line.startswith(">")]
    report = pd.DataFrame(headers, columns=["qseqid"])
    cog_report = pd.read_csv(f'{out}/Annotation/{sample}/COG_report.tsv', sep='\t', low_memory=False)
    cog_report = cog_report.groupby('qseqid')[cog_report.columns.tolist()[1:]].first().reset_index()
    report = pd.merge(report, cog_report, on='qseqid', how='left')
    report = report.groupby('qseqid')[report.columns.tolist()[1:]].first().reset_index()
    report = report[report['DB ID'].str.startswith('COG') == True].rename(columns={'DB ID': 'COG ID'})
    report = pd.merge(
        pd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t', low_memory=False),
        report, on='qseqid', how='outer')
    rename_cols = blast_cols + ['EC number']
    report = report.rename(columns={
        **{f'{col}_x': f'{col} (UPIMAPI)' for col in rename_cols},
        **{f'{col}_y': f'{col} (reCOGnizer)' for col in rename_cols}})
    report['Contig'] = report['qseqid'].apply(lambda x: x.split('_')[1])
    mg_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'dna')]['Name'].tolist()
    mt_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'mrna')]['Name'].tolist()
    mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()
    if len(mg_names) > 0:
        mg_counts = pd.read_csv(f'{out}/Quantification/{sample}_mg_norm.tsv', sep='\t')
        mg_counts['Contig'] = mg_counts['Contig'].apply(lambda x: x.split('_')[1])
        report = pd.merge(report, mg_counts, on='Contig', how='left')
        mg_preport = pd.merge(mg_preport, report[['Entry'] + mg_names], on='Entry', how='outer')
    if len(mt_names) > 0:
        report = pd.merge(
            report, pd.read_csv(f'{out}/Quantification/{sample}_mt_norm.tsv', sep='\t', names=[
                'qseqid'] + mt_names), on='qseqid', how='left')
        mt_preport = pd.merge(mt_preport, report[['Entry'] + mt_names], on='Entry', how='outer')
        # MT readcounts come normalized by gene size, and therefore not fit for DE (as they are floats).
        # DE input are the non-normalized readcounts.
        de_input = pd.merge(
            de_input, pd.read_csv(f'{out}/Quantification/{sample}_mt.readcounts', sep='\t').rename(
                columns={'Gene': 'Entry'}), on='Entry', how='outer')
    if len(mp_names) > 0:
        spectracounts = pd.read_csv(f'{out}/Metaproteomics/{sample}/spectracounts.tsv', sep='\t')
        spectracounts.rename(columns={'Main Accession': 'qseqid'}, inplace=True)
        report = pd.merge(report, spectracounts, on='qseqid', how='left')
        mp_preport = pd.merge(mp_preport, report[['Entry'] + mp_names], on='Entry', how='outer')
    report[mg_names + mt_names + mp_names] = report[mg_names + mt_names + mp_names].fillna(
        value=0).astype(float).astype(int)
    report.to_csv(f'{out}/MOSCA_{sample}_Protein_Report.tsv', sep='\t', index=False)
    return report, mg_preport, mt_preport, mp_preport, de_input


def make_protein_reports(out, exps, max_lines=1000000):
    mg_report = mt_report = mp_report = de_input = pd.DataFrame(columns=['Entry'])
    writer = pd.ExcelWriter(f'{out}/MOSCA_Protein_Report.xlsx', engine='xlsxwriter')
    for sample in set(exps['Sample']):
        report, mg_report, mt_report, mp_report, de_input = make_protein_report(
            out, exps, sample, mg_report, mt_report, mp_report, de_input)
        timed_message(f'Writing Protein Report for sample: {sample}.')
        if len(report) < max_lines:
            report.to_excel(writer, sheet_name=sample, index=False)
        else:
            k = 1
            for i in range(0, len(report), max_lines):
                j = min(i + max_lines, len(report))
                report.iloc[i:(i + j)].to_excel(writer, sheet_name=f'{sample} ({k})', index=False)
                k += 1
    writer.close()
    # Write quantification matrices to normalize all together
    timed_message('Writing quantification matrices.')
    if len(mg_report) > 0:
        mg_report[mg_report.columns.tolist()[1:]] = mg_report[mg_report.columns.tolist()[1:]].astype(float)
        mg_report = mg_report.groupby('Entry')[mg_report.columns.tolist()[1:]].sum().reset_index()
        mg_report.to_csv(f'{out}/Quantification/mg_entry_quant.tsv', sep='\t', index=False)
    if len(mt_report) > 0:
        mt_report[mt_report.columns.tolist()[1:]] = mt_report[mt_report.columns.tolist()[1:]].astype(float)
        mt_report = mt_report.groupby('Entry')[mt_report.columns.tolist()[1:]].sum().reset_index()
        mt_report.to_csv(f'{out}/Quantification/mt_entry_quant.tsv', sep='\t', index=False)
        # MT readcounts come normalized by gene size, and therefore not fit for DE (as they are floats).
        # DE input are the non-normalized readcounts.
        de_input.to_csv(f'{out}/Quantification/dea_input.tsv', sep='\t', index=False)
    if len(mp_report) > 0:
        mp_report[mp_report.columns.tolist()[1:]] = mp_report[mp_report.columns.tolist()[1:]].astype(float)
        mp_report = mp_report.groupby('Entry')[mp_report.columns.tolist()[1:]].sum().reset_index()
        mp_report[mp_report[mp_report.columns.tolist()[1:]].isnull().sum(
            axis=1) < len(mp_report.columns.tolist()[1:])].drop_duplicates().to_csv(
            f'{out}/Metaproteomics/mp_entry_quant', sep='\t', index=False)
        # in the case of MP, there is no normalization by protein size. So it can be done this way
        shutil.copyfile(f'{out}/Metaproteomics/mp_entry_quant', f'{out}/Quantification/dea_input.tsv')


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')
    make_protein_reports(snakemake.params.output, exps)


if __name__ == '__main__':
    run()
