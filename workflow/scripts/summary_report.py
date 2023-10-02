# -*- coding: utf-8 -*-
"""
General report construction and export

By João Sequeira

Oct 2019
"""

from glob import glob
import pandas as pd
from zipfile import ZipFile
from mosca_tools import run_pipe_command, parse_fastqc_report, count_on_file, timed_message


class Reporter:
    def __init__(self):
        self.report = pd.DataFrame(columns=[
            'Initial reads', 'Qual trim params', 'Final reads', '# contigs', 'N50', 'Reads aligned (%)',
            '# high-quality MAGs', '# medium-quality MAGs', '# low-quality MAGs', '# genes',
            '# annotations (UPIMAPI)', '# annotations (reCOGnizer)', '# differentially expressed'])

    def write_technical_report(self, output):
        """
        Writes the report with the softwares used by MOSCA and respective versions to a file
        param: output: str - path to output file
        """
        timed_message('Writting technical report.')
        conda_list = run_pipe_command('conda list', output='PIPE', print_message=False).split('\n')[2:]
        lines = [line.split() for line in conda_list]
        lines[0] = lines[0][1:]
        df = pd.DataFrame(lines, columns=lines.pop(0)).set_index('Name')
        df[['Version']].to_csv(output, sep='\t')

    def info_from_preprocessing(self, out_dir):
        timed_message('Obtaining info from preprocessing.')
        reports = glob(f'{out_dir}/Preprocess/FastQC/*/fastqc_data.txt')
        for file in reports:
            timed_message(f'Obtaining info from: {file}')
            if 'noadapters' in file or 'norrna' in file:
                continue
            name = file.split('/')[-2].split('_R')[0].split('_forward')[0].split('_reverse')[0].split('_trimmed_')[-1]
            self.report = pd.concat([self.report, pd.Series(name=name, dtype='object')])
            if 'quality_trimmed' in file:
                self.report.loc[name, 'Final reads'] = parse_fastqc_report(file)['Basic Statistics'][1].loc[
                    'Total Sequences', 'Value']
            else:
                self.report.loc[name, 'Initial reads'] = parse_fastqc_report(file)['Basic Statistics'][1].loc[
                    'Total Sequences', 'Value']
            with open(f'{out_dir}/Preprocess/Trimmomatic/{name}_quality_params.txt') as f:
                self.report.loc[name, 'Qual trim params'] = ';'.join(f.read().split('\n'))

    def info_from_assembly(self, out_dir):
        timed_message('Obtaining info from assembly.')
        reports = glob(f'{out_dir}/Assembly/*/quality_control/report.tsv')
        for file in reports:
            timed_message(f'Obtaining info from: {file}')
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
            timed_message(f'Obtaining info from: {file}')
            sample = file.split('/')[-2]
            if sample not in self.report.index:
                self.report = pd.concat([self.report, pd.Series(name=sample, dtype='object')])
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
            timed_message(f'Obtaining info from: {file}')
            sample = file.split('/')[-2]
            if sample not in self.report.index:
                self.report = pd.concat([self.report, pd.Series(name=sample, dtype='object')])
            self.report['# genes'] = count_on_file('>', file)
        for upi_res in upimapi_res:
            timed_message(f'Obtaining info from: {upi_res}')
            sample = file.split('/')[-2]
            self.report.loc[sample, '# annotations (UPIMAPI)'] = pd.read_csv(
                upi_res, sep='\t', low_memory=False)['qseqid'].unique().sum()
        for recog_res in recognizer_res:
            timed_message(f'Obtaining info from: {recog_res}')
            sample = file.split('/')[-2]
            self.report.loc[sample, '# annotations (reCOGnizer)'] = pd.read_csv(
                recog_res, sep='\t', low_memory=False)['qseqid'].unique().sum()

    def info_from_mt_quantification(self, out_dir):
        timed_message('Obtaining info from mt quantification.')
        reports = glob(f'{out_dir}/Quantification/*.log')
        for file in reports:
            timed_message(f'Obtaining info from: {file}')
            name = file.split('/')[-1].split('.log')[0]
            if name not in self.report.index:
                self.report = pd.concat([self.report, pd.Series(name=name, dtype='object')])
            with open(file) as f:
                self.report.loc[name, 'Reads aligned (%)'] = f.readlines()[-1].split('%')[0]

    def info_from_differential_expression(self, out_dir, cutoff=0.01):
        timed_message('Obtaining info from differential expression.')
        reports = glob(f'{out_dir}/Quantification/*/condition_treated_results.tsv')
        for file in reports:
            timed_message(f'Obtaining info from: {file}')
            sample = file.split('/')[-2]
            if sample not in self.report.index:
                self.report = pd.concat([self.report, pd.Series(name=sample, dtype='object')])
            de_results = pd.read_csv(f'{out_dir}/Quantification/{sample}/condition_treated_results.tsv', sep='\t')
            self.report.loc[sample, '# differentially expressed'] = (de_results['padj'] < cutoff).sum()

    def zip_outputs(self, out_dir):
        files_n_folders = {
            'fastqc_reports': [file for file in glob(f'{out_dir}/Preprocess/FastQC/*.html') if (
                'noadapters' not in file and 'norrna' not in file)],
            'assembly_reports': glob(f'{out_dir}/Assembly/*/quality_control/report.tsv'),
            'taxonomy_kronas': glob(f'{out_dir}/*_tax.html'),
            'functional_kronas': glob(f'{out_dir}/*_fun.html'),
            'de_plots': glob(f'{out_dir}/Quantification/*/*.jpeg'),
            'kegg_maps': glob(f'{out_dir}/KEGG_maps/*.png'),
            'main_reports': [f'{out_dir}/{filename}' for filename in ['MOSCA_Protein_Report.xlsx',
                    'MOSCA_Entry_Report.xlsx', 'MOSCA_General_Report.tsv', 'technical_report.tsv']]}
        with ZipFile(f'{out_dir}/MOSCA_results.zip', 'w') as archive:
            for k, v in files_n_folders.items():
                for file in v:
                    prefix = file.split('/')[-3] + '_' if k == 'assembly_reports' else (
                        file.split('/')[-2] + '_' if k in ['de_plots', 'kegg_maps'] else '')
                    archive.write(file, arcname=f'{k}/{prefix}{file.split("/")[-1]}')

    def run(self):
        timed_message('Writting final reports.')
        self.write_technical_report(f'{snakemake.params.output}/technical_report.tsv')
        self.info_from_preprocessing(snakemake.params.output)
        self.info_from_assembly(snakemake.params.output)
        self.info_from_binning(snakemake.params.output)
        self.info_from_annotation(snakemake.params.output)
        self.info_from_mt_quantification(snakemake.params.output)
        self.info_from_differential_expression(snakemake.params.output)
        self.report.dropna(how='all', axis=1).to_csv(f'{snakemake.params.output}/MOSCA_General_Report.tsv', sep='\t')
        self.zip_outputs(snakemake.params.output)


if __name__ == '__main__':
    Reporter().run()
