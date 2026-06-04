# -*- coding: utf-8 -*-
"""
Summary report construction and export

By João Sequeira

Oct 2019
"""

from mosca_tools import run_pipe_command, parse_fastqc_report, count_on_file, timed_message
from io import StringIO
from yaml import safe_load 
import pandas as pd
from zipfile import ZipFile
from glob import glob


def get_env_info(env_path=''):
    conda_list = run_pipe_command(f"conda list{' -p ' + env_path if env_path else ''}", output='PIPE', print_message=False)
    df = pd.read_csv(StringIO(conda_list), skiprows=3, sep='\s+', header=None, names=['Name', 'Version', 'Build', 'Channel'])
    if env_path:
        with open(f'{env_path}.yaml') as file_data:
            yaml_data = safe_load(file_data)
            return yaml_data['name'], df
    return 'base', df
    

def write_versions_report(output):
    """
    Writes the report with the softwares used by MOSCA and respective versions to a file
    param: output: str - path to output file
    """
    timed_message('Writting technical report.')
    envs = [f for f in glob('.snakemake/conda/*') if not ('.' in f)]

    writer = pd.ExcelWriter(output, engine='xlsxwriter')

    for env in [''] + envs:         # '' is for the base environment
        name, df = get_env_info(env)
        df.to_excel(writer, sheet_name=name if env else 'base', index=False)

    writer.close()


def add_preprocessing(report, sample2name, out_dir):
    timed_message('Processing preprocessing.')
    reports = glob(f'{out_dir}/Preprocess/FastQC/*/fastqc_data.txt')

    for file in reports:
        if 'noadapters' in file or 'norrna' in file:
            continue

        name = file.split('/')[-2].split('_R')[0].split('_forward')[0].split('_reverse')[0].split('_trimmed_')[-1]

        stats = parse_fastqc_report(file)['Basic Statistics'][1]
        val = stats.loc['Total Sequences', 'Value']

        if 'quality_trimmed' in file:
            report.loc[name, 'Final reads'] = val
        else:
            report.loc[name, 'Initial reads'] = val

        with open(f'{out_dir}/Preprocess/Trimmomatic/{name}_quality_params.txt') as f:
            report.loc[name, 'Qual trim params'] = ';'.join([x for x in f.read().split('\n') if x])

    return report, sample2name


def add_assembly(report, sample2name, out_dir):
    timed_message('Processing assembly.')
    reports = glob(f'{out_dir}/Assembly/*/quality_control/report.tsv')

    for file in reports:
        sample = file.split('/')[-3]

        report = pd.concat([report, pd.Series(name=sample, dtype='object')])

        data = pd.read_csv(file, sep='\t', index_col='Assembly')

        name = sample2name.get(sample, sample)

        report.loc[name, ['# contigs', 'N50', 'Reads aligned (%)']] = (
            int(data.loc['# contigs', 'contigs']),
            int(data.loc['N50', 'contigs']),
            data.loc['Reads aligned (%)', 'contigs']
        )

    return report, sample2name


def add_binning(report, sample2name, out_dir):
    timed_message('Processing binning.')
    reports = glob(f'{out_dir}/Binning/*/checkm.tsv')

    for file in reports:
        sample = file.split('/')[-2]
        data = pd.read_csv(file, sep='\t')

        name = sample2name.get(sample, sample)

        report.loc[name, [
            '# high-qual MAGs',
            '# medium-qual MAGs',
            '# low-qual MAGs'
        ]] = (
            ((data['Completeness'] >= 90) & (data['Contamination'] <= 5)).sum(),
            ((data['Completeness'] >= 50) & (data['Completeness'] < 90) & (data['Contamination'] <= 10)).sum(),
            ((data['Completeness'] < 50) & (data['Contamination'] <= 10)).sum()
        )

    return report


def add_annotation(report, sample2name, out_dir):
    timed_message('Processing annotation.')

    fastas = glob(f'{out_dir}/Annotation/*/fgs.faa')
    upimapi = glob(f'{out_dir}/Annotation/*/UPIMAPI_results.tsv')
    recognizer = glob(f'{out_dir}/Annotation/*/reCOGnizer_results.tsv')

    for file in fastas:
        sample = file.split('/')[-2]
        name = sample2name.get(sample, sample)
        report.loc[name, '# genes'] = count_on_file('>', file)

    for file in upimapi:
        sample = file.split('/')[-2]
        name = sample2name.get(sample, sample)
        report.loc[name, '# annotations (UPIMAPI)'] = len(pd.read_csv(file, sep='\t', low_memory=False)['qseqid'].unique())

    for file in recognizer:
        sample = file.split('/')[-2]
        name = sample2name.get(sample, sample)
        report.loc[name, '# annotations (reCOGnizer)'] = len(pd.read_csv(file, sep='\t', low_memory=False)['qseqid'].unique())

    return report


def add_quantification(report, sample2name, out_dir):
    timed_message('Processing quantification.')

    logs = glob(f'{out_dir}/Quantification/*.log')

    for file in logs:
        name = file.split('/')[-1].split('.log')[0]
        with open(file) as f:
            report.loc[name, 'Reads aligned (%)'] = f.readlines()[-1].split('%')[0]

    return report


def add_de(report, sample2name, out_dir, cutoff=0.01, mp=False):
    timed_message('Processing DE analysis.')

    file = f'{out_dir}/DE_analysis/condition_treated_results.tsv'
    de = pd.read_csv(file, sep='\t', index_col=0)

    report['# differentially expressed'] = (
        (de['pvalue'] < cutoff) &
        (de['FDR' if mp else 'padj'] < cutoff)
    ).sum()

    return report


def zip_outputs(out_dir):
    timed_message(f'Zipping results: {out_dir}/MOSCA_results.zip')

    files = {
        'fastqc_reports': glob(f'{out_dir}/Preprocess/FastQC/*.html'),
        'assembly_reports': glob(f'{out_dir}/Assembly/*/quality_control/report.tsv'),
        'taxonomy_kronas': glob(f'{out_dir}/kronas/*_tax.html'),
        'functional_kronas': glob(f'{out_dir}/kronas/*_fun.html'),
        'de_plots': glob(f'{out_dir}/DE_analysis/*.jpeg'),
        'kegg_maps': glob(f'{out_dir}/KEGG_maps/*.png'),
    }

    with ZipFile(f'{out_dir}/MOSCA_results.zip', 'w') as archive:
        for k, v in files.items():
            for f in v:
                archive.write(f, arcname=f'{k}/{f.split("/")[-1]}')

def run():
    out = snakemake.params.output

    report, sample2name = pd.DataFrame(), {}

    write_versions_report(f'{out}/MOSCA_Versions_Report.xlsx')

    exps = pd.read_csv(f'{out}/exps.tsv', sep='\t')

    for sample in exps['Sample'].unique():
        sample2name[sample] = exps[exps['Sample'] == sample]['Name'].tolist()

    report, sample2name = add_preprocessing(report, sample2name, out)
    report, sample2name = add_assembly(report, sample2name, out)
    report = add_binning(report, sample2name, out)
    report = add_annotation(report, sample2name, out)
    report = add_quantification(report, sample2name, out)

    if snakemake.params.has_expression_data:
        report = add_de(
            report, sample2name, out,
            cutoff=snakemake.params.cutoff,
            mp='protein' in exps['Data type'].tolist()
        )

    cols = [
        'Initial reads', 'Qual trim params', 'Final reads',
        '# genes', '# annotations (UPIMAPI)',
        '# annotations (reCOGnizer)', 'Reads aligned (%)'
    ]

    if snakemake.params.has_expression_data:
        cols.append('# differentially expressed')

    report[cols].to_csv(f'{out}/MOSCA_Summary_Report.tsv', sep='\t')

    zip_outputs(out)


if __name__ == '__main__':
    run()