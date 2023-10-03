# -*- coding: utf-8 -*-
"""
General report construction and export

By João Sequeira

Oct 2019
"""

from glob import glob
import pandas as pd
from zipfile import ZipFile

import yaml

from mosca_tools import run_pipe_command, parse_fastqc_report, count_on_file, timed_message


def get_env_info(yaml_file):
    hash_string = yaml_file.split('/')[-1].split('.')[0]
    with open(yaml_file) as file_data:
        yaml_data = yaml.safe_load(file_data)
    deps = [dep.split('=') for dep in yaml_data['dependencies'] if type(dep) == str]
    pip_deps = [dep for dep in yaml_data['dependencies'] if type(dep) == dict]
    pip_deps = [dep.split('==') for dep in (pip_deps[0]['pip'] if len(pip_deps) > 0 else [])]
    versions = pd.DataFrame.from_dict(
        {dep[0]: dep[1:] for dep in deps + pip_deps}, orient='index', columns=['Version', 'Build'])
    versions.index.set_names('Name', inplace=True)
    return hash_string, yaml_data['name'], yaml_data['channels'], versions


def write_versions_report(output):
    """
    Writes the report with the softwares used by MOSCA and respective versions to a file
    param: output: str - path to output file
    """
    timed_message('Writting technical report.')
    yamls = glob('.snakemake/conda/*.yaml')
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    conda_list = run_pipe_command('conda list', output='PIPE', print_message=False)
    values = [line.split() for line in conda_list.split('\n')[2:]][:-1]
    df = pd.DataFrame(values[1:], columns=values[0][1:])
    df.to_excel(writer, sheet_name='Base env', index=False)
    for file in yamls:
        hash_string, name, channels, versions = get_env_info(file)
        versions.to_excel(writer, sheet_name=name)
    writer.close()


class Reporter:
    def __init__(self):
        self.report = pd.DataFrame()


    def info_from_preprocessing(self, out_dir):
        timed_message('Obtaining info from preprocessing.')
        reports = glob(f'{out_dir}/Preprocess/FastQC/*/fastqc_data.txt')
        for file in reports:
            if 'noadapters' in file or 'norrna' in file:
                continue
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-3:])}')
            name = file.split('/')[-2].split('_R')[0].split('_forward')[0].split('_reverse')[0].split('_trimmed_')[-1]
            if 'quality_trimmed' in file:
                self.report.loc[name, 'Final reads'] = parse_fastqc_report(file)['Basic Statistics'][1].loc[
                    'Total Sequences', 'Value']
            else:
                self.report.loc[name, 'Initial reads'] = parse_fastqc_report(file)['Basic Statistics'][1].loc[
                    'Total Sequences', 'Value']
            with open(f'{out_dir}/Preprocess/Trimmomatic/{name}_quality_params.txt') as f:
                self.report.loc[name, 'Qual trim params'] = ';'.join([x for x in f.read().split('\n') if len(x) > 0])

    def info_from_assembly(self, out_dir):
        timed_message('Obtaining info from assembly.')
        reports = glob(f'{out_dir}/Assembly/*/quality_control/report.tsv')
        for file in reports:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-4:])}')
            sample = file.split('/')[-3]
            self.report = pd.concat([self.report, pd.Series(name=sample, dtype='object')])
            data = pd.read_csv(file, sep='\t', index_col='Assembly')
            self.report.loc[sample, ['# contigs', 'N50', 'Reads aligned (%)']] = (
                int(data.loc['# contigs', 'contigs']),
                int(data.loc['N50', 'contigs']),
                data.loc['Reads aligned (%)', 'contigs'])

    def info_from_binning(self, out_dir):
        timed_message('Obtaining info from binning.')
        reports = glob(f'{out_dir}/Binning/*/checkm.tsv')
        for file in reports:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-3:])}')
            sample = file.split('/')[-2]
            data = pd.read_csv(file, sep='\t')
            self.report.loc[sample, ['# high-qual MAGs', '# medium-qual MAGs', '# low-qual MAGs']] = (
                ((data['Completeness'] >= 90) & (data['Contamination'] <= 5)).sum(),
                ((data['Completeness'] >= 50) & (data['Completeness'] < 90) & (data['Contamination'] <= 10)).sum(),
                ((data['Completeness'] < 50) & (data['Contamination'] <= 10)).sum())

    def info_from_annotation(self, out_dir):
        timed_message('Obtaining info from annotation.')
        fastas = glob(f'{out_dir}/Annotation/*/fgs.faa')
        upimapi_res = glob(f'{out_dir}/Annotation/*/UPIMAPI_results.tsv')
        recognizer_res = glob(f'{out_dir}/Annotation/*/reCOGnizer_results.tsv')
        for file in fastas:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-3:])}')
            sample = file.split('/')[-2]
            self.report.loc[sample, '# genes'] = count_on_file('>', file)
        for file in upimapi_res:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-3:])}')
            sample = file.split('/')[-2]
            self.report.loc[sample, '# annotations (UPIMAPI)'] = pd.read_csv(
                file, sep='\t', low_memory=False)['qseqid'].unique().sum()
        for file in recognizer_res:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-3:])}')
            sample = file.split('/')[-2]
            self.report.loc[sample, '# annotations (reCOGnizer)'] = pd.read_csv(
                file, sep='\t', low_memory=False)['qseqid'].unique().sum()

    def info_from_quantification(self, out_dir):
        timed_message('Obtaining info from mt quantification.')
        reports = glob(f'{out_dir}/Quantification/*.log')
        for file in reports:
            timed_message(f'Obtaining info from: {"/".join(file.split("/")[-2:])}')
            name = file.split('/')[-1].split('.log')[0]
            with open(file) as f:
                self.report.loc[name, 'Reads aligned (%)'] = f.readlines()[-1].split('%')[0]

    def info_from_differential_expression(self, out_dir, cutoff=0.01, mp=False):
        file = f'{out_dir}/DE_analysis/condition_treated_results.tsv'
        timed_message(f'Obtaining info from: {"/".join(file.split("/")[-2:])}')
        de_results = pd.read_csv(file, sep='\t', index_col=0)
        self.report['# differentially expressed'] = (       # cutoff applies to both pvalue and FDR/padj
                (de_results['pvalue'] < cutoff) & (de_results['FDR' if mp else 'padj'] < cutoff)).sum()

    def zip_outputs(self, out_dir):
        files_n_folders = {
            'fastqc_reports': [file for file in glob(f'{out_dir}/Preprocess/FastQC/*.html') if (
                'noadapters' not in file and 'norrna' not in file)],
            'assembly_reports': glob(f'{out_dir}/Assembly/*/quality_control/report.tsv'),
            'taxonomy_kronas': glob(f'{out_dir}/kronas/*_tax.html'),
            'functional_kronas': glob(f'{out_dir}/kronas/*_fun.html'),
            'de_plots': glob(f'{out_dir}/DE_analysis/*.jpeg'),
            'kegg_maps': glob(f'{out_dir}/KEGG_maps/*.png'),
            'main_reports': [f'{out_dir}/{filename}' for filename in [
                'MOSCA_Protein_Report.xlsx', 'MOSCA_Entry_Report.xlsx', 'MOSCA_General_Report.tsv',
                'technical_report.tsv']]}
        with ZipFile(f'{out_dir}/MOSCA_results.zip', 'w') as archive:
            for k, v in files_n_folders.items():
                for file in v:
                    prefix = file.split('/')[-3] + '_' if k == 'assembly_reports' else ''   # get sample name to distinguish assembly reports
                    archive.write(file, arcname=f'{k}/{prefix}{file.split("/")[-1]}')

    def run(self):
        timed_message('Writting final reports.')
        write_versions_report(f'{snakemake.params.output}/MOSCA_Versions_Report.xlsx')
        self.info_from_preprocessing(snakemake.params.output)
        self.info_from_assembly(snakemake.params.output)
        self.info_from_binning(snakemake.params.output)
        self.info_from_annotation(snakemake.params.output)
        self.info_from_quantification(snakemake.params.output)
        exps = pd.read_csv(f'{snakemake.params.output}/exps.tsv', sep='\t')
        self.info_from_differential_expression(
            snakemake.params.output, cutoff=snakemake.params.cutoff, mp='protein' in exps['Data type'].tolist())
        self.report.to_csv(f'{snakemake.params.output}/MOSCA_General_Report.tsv', sep='\t')
        self.zip_outputs(snakemake.params.output)


if __name__ == '__main__':
    Reporter().run()
