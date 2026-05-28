"""
MOSCA's script for producing General report

By João Sequeira

Dec 2022
"""

from mosca_tools import timed_message, blast_cols
import pandas as pd
import numpy as np
import dask.dataframe as dd
import shutil
from concurrent.futures import ThreadPoolExecutor

functional_columns = [
    'General functional category', 'Functional category', 'Protein description', 'COG ID', 'EC number (reCOGnizer)']

def process_quantification(data_type, did_assembly, names, out, sample):
    """
    Read counts (reads or spectra) depending on the type of data (mg, mt, mp).
    """
    if data_type == 'mg':
        filepath = f'{out}/Quantification/{sample}_mg_norm.tsv' if did_assembly else f'{out}/Quantification/{sample}_mg.readcounts'
        counts = pd.read_csv(filepath, sep='\t')            # TODO - check if when no assembly the first column is named qseqid in the norm file
        if did_assembly:
            counts['Contig'] = counts['Contig'].str.split('_').str[1]
        return counts.set_index('Contig' if did_assembly else 'qseqid')
    if data_type == 'mt':
        norm_filepath = f'{out}/Quantification/{sample}_mt_norm.tsv' if did_assembly else f'{out}/Quantification/{sample}_mt.readcounts'
        counts = pd.read_csv(norm_filepath, sep='\t')
        return counts.rename(columns={'Gene': 'qseqid'}).set_index('qseqid')
    if data_type == 'mp':
        counts = pd.read_csv(f'{out}/Metaproteomics/{sample}_mp.spectracounts', sep='\t')
        return counts.rename(columns={'Main Accession': 'qseqid'}).set_index('qseqid')
    return ValueError(f'Unknown data type: {data_type}')


def make_general_report(out, exps, sample, mg_preport, mt_preport, mp_preport, de_input, did_assembly=True):
    timed_message(f'Joining data for sample: {sample}.')
    print('Reading gene calling headers.')
    with open(f'{out}/Annotation/{sample}/fgs.faa') as f:
        headers = [line.strip()[1:] for line in f if line.startswith(">")]

    report = pd.DataFrame(headers, columns=["qseqid"]).set_index('qseqid')

    print('Reading reCOGnizer results.')
    cog_report = dd.read_csv(f'{out}/Annotation/{sample}/COG_report.tsv', sep='\t', dtype=str)
    cog_report = cog_report[cog_report['DB ID'].str.startswith('COG') == True].rename(columns={'DB ID': 'COG ID'})
    cog_report = cog_report.groupby('qseqid').first().compute()
    report = pd.merge(report, cog_report, left_index=True, right_index=True, how='left')

    print('Reading UPIMAPI results.')
    upimapi_results = dd.read_csv(f'{out}/Annotation/{sample}/UPIMAPI_results.tsv', sep='\t', dtype=str).compute().set_index('qseqid')
    report = pd.merge(upimapi_results, report, left_index=True, right_index=True, how='outer')

    print('Formatting names of columns.')
    rename_cols = blast_cols + ['EC number']
    report = report.rename(columns={**{f'{col}_x': f'{col} (UPIMAPI)' for col in rename_cols},
                                    **{f'{col}_y': f'{col} (reCOGnizer)' for col in rename_cols}})

    if did_assembly:
        report['Contig'] = report.index.to_series().str.split('_').str[1]

    mg_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'dna')]['Name'].tolist()
    mt_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'mrna')]['Name'].tolist()
    mp_names = exps[(exps['Sample'] == sample) & (exps['Data type'] == 'protein')]['Name'].tolist()

    with ThreadPoolExecutor() as executor:
        futures = {}
        if mg_names:
            futures['mg'] = executor.submit(process_quantification, 'mg', did_assembly, mg_names, out, sample)
        if mt_names:
            futures['mt'] = executor.submit(process_quantification, 'mt', did_assembly, mt_names, out, sample)
        if mp_names:
            futures['mp'] = executor.submit(process_quantification, 'mp', did_assembly, mp_names, out, sample)
        results = {key: future.result() for key, future in futures.items()}

    for dtype, result in results.items():
        if dtype == 'mg' and did_assembly:
            report = pd.merge(report, result, left_on='Contig', right_index=True, how='left')
        else:
            report = pd.merge(report, result, left_index=True, right_index=True, how='left')
    
    dfs=[]; 
    for dtype, result in results.items():
        result.index.name='Entry'
        if dtype=='mg':
            mg_preport=pd.merge(mg_preport,result,on='Entry',how='outer')
        elif dtype=='mt':
            mt_preport=pd.merge(mt_preport,result,on='Entry',how='outer')
            for name in mt_names:
                dfs.append(pd.read_csv(f'{out}/Quantification/{name}.readcounts', sep='\t', names=['Entry',name]).set_index('Entry'))
        elif dtype=='mp':
            mp_preport=pd.merge(mp_preport,result,on='Entry',how='outer')

    if dfs:
        de_input=pd.concat(dfs,axis=1,join='outer').reset_index()
        de_input['Entry']=de_input['Entry'].map(upimapi_results['sseqid']).fillna(de_input['Entry'])

    report[mg_names + mt_names + mp_names] = report[mg_names + mt_names + mp_names].fillna(0).astype(float)     # astype(float).astype(int) avoids "ValueError: invalid literal for int() with base 10: '2.0'"
    report.to_csv(f'{out}/MOSCA_{sample}_General_Report.tsv', sep='\t', index=False)
    return report, mg_preport, mt_preport, mp_preport, de_input


def make_general_reports(out, exps, max_lines=1000000, did_assembly=True):
    mg_report = mt_report = mp_report = pd.DataFrame(columns=['Entry'])
    de_input = pd.DataFrame()
    writer = pd.ExcelWriter(f'{out}/MOSCA_General_Report.xlsx', engine='xlsxwriter')

    for sample in set(exps['Sample']):
        report, mg_report, mt_report, mp_report, de_input = make_general_report(
            out, exps, sample, mg_report, mt_report, mp_report, de_input, did_assembly=did_assembly)
        timed_message(f'Writing General Report for sample: {sample}')
        if len(report) < max_lines:
            report.to_excel(writer, sheet_name=sample, index=False)
        else:
            for k, chunk in enumerate(np.array_split(report, len(report) // max_lines)):
                chunk.to_excel(writer, sheet_name=f'{sample} ({k + 1})', index=False)
    writer.close()

    if len(de_input) > 0:
        de_input = de_input.groupby('Entry').sum().reset_index()
        de_input[de_input.columns.tolist()[1:]] = de_input[de_input.columns.tolist()[1:]].fillna(0).astype(int)
        de_input.to_csv(f'{out}/Quantification/dea_input.tsv', sep='\t', index=False)

    timed_message('Writing quantification matrices.')
    if 'dna' in exps['Data type'].values:
        mg_report.iloc[:, 1:] = mg_report.iloc[:, 1:].astype(float)
        mg_report = mg_report.groupby('Entry').sum()
        mg_report.to_csv(f'{out}/Quantification/mg_entry_quant.tsv', sep='\t')
    if 'mrna' in exps['Data type'].values:
        mt_report.iloc[:, 1:] = mt_report.iloc[:, 1:].astype(float)
        mt_report = mt_report.groupby('Entry').sum()
        mt_report.to_csv(f'{out}/Quantification/mt_entry_quant.tsv', sep='\t')
    if 'protein' in exps['Data type'].values:
        mp_report.iloc[:, 1:] = mp_report.iloc[:, 1:].astype(float)
        mp_report = mp_report.groupby('Entry').sum()
        mp_report = mp_report.drop_duplicates().dropna(subset=mp_report.columns[1:])
        mp_report.to_csv(f'{out}/Metaproteomics/mp_entry_quant.tsv', sep='\t')


def run():
    exps = pd.read_csv(snakemake.params.exps, sep='\t')
    make_general_reports(snakemake.params.output, exps, did_assembly=snakemake.params.did_assembly)


if __name__ == '__main__':
    run()
