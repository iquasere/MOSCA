# -*- coding: utf-8 -*-
'''
General report construction and export

By João Sequeira

Oct 2019
'''

import argparse
import glob
import pandas as pd
from mosca_tools import parse_blast, run_pipe_command, parse_fastqc_report, count_on_file
from zipfile import ZipFile


class Reporter:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA's technical and quality control reports")
        parser.add_argument("-e", "--experiments", type=str, required=True,
                            help="Experiments file")
        parser.add_argument("-o", "--output", type=str, help="Output directory"),
        parser.add_argument("-ldir", "--lists-directory", type=str, help="Directory with lists for Reporter")
        parser.add_argument("-if", "--input-format", type=str, default='tsv', choices=['tsv', 'excel'])
        args = parser.parse_args()
        args.output = args.output.rstrip('/')
        return args

    '''
    Input:
        output: str - filename to write tools and respective versions
    Output:
        a file named [output] will be written with information concerning
        the softwares used by MOSCA and respective versions
    '''

    def write_technical_report(self, output):  # TODO - add proteomics software that cannot be installed with conda
        conda_list = run_pipe_command('conda list', output='PIPE').split('\n')[2:]
        lines = [line.split() for line in conda_list]
        lines[0] = lines[0][1:]
        df = pd.DataFrame(lines, columns=lines.pop(0)).set_index('Name')
        df[['Version']].to_csv(output, sep='\t')

    def initialize_report(self, reporter_columns):
        print('Initializing Report')
        self.report = pd.DataFrame(columns=reporter_columns)

    '''
    Input: 
        name: str - name of MG/MT sample
        prefix: str - [Initial quality assessment], [Before quality trimming] or
        [After quality trimming] - the step that concerns the FastQC report
        performed_rrna_removal: bool - True if rRNA removal was performed, False
        otherwise. Should be False for MG, and True for MT
    Output:
        self.report will be updated with information from the report, on the line
        named 'name', and the columns that start with 'prefix'
    '''

    def info_from_fastqc(self, output_dir, name, prefix, prefix2terms, fastq_columns, initial_file=None):
        if initial_file is not None:
            if '_R' in initial_file:
                original_name = initial_file.split('/')[-1].split('_R')[0]
            else:
                original_name = initial_file.split('/')[-1].split('.f')[0]
        else:
            original_name = name
        reports = [parse_fastqc_report('{}/Preprocess/FastQC/{}{}_{}_fastqc/fastqc_data.txt'.format(
            output_dir, prefix2terms[prefix][0], original_name, prefix2terms[prefix][i])) for i in [1, 2]]
        self.report.loc[name, '{} # of reads'.format(prefix)] = reports[0]['Basic Statistics'][1].loc[
            'Total Sequences']['Value']

        for column in fastq_columns:
            for i in range(2):
                if column not in reports[i].keys():
                    reports[i][column] = (
                    'Not available', None)  # only the not available matters. And nothing else matters!...
            if reports[0][column][0] == reports[1][column][0]:
                self.report.loc[name, '{} {}'.format(prefix, column)] = reports[0][column][0]
            else:
                self.report.loc[name, '{} {}'.format(prefix, column)] = (
                    '{} (forward) {} (reverse)'.format(
                        reports[0][column][0], reports[1][column][0]))

    def info_from_preprocessing(self, output_dir, name, input_file, fastq_columns, performed_rrna_removal=False):
        print('Retrieving preprocessing information for dataset: ' + name)
        if name not in self.report.index:
            self.report = self.report.append(pd.Series(name=name))
        self.report.loc[name] = self.report.loc[name].fillna(value='')

        adapter_files = open('{}/Preprocess/Trimmomatic/{}_adapters.txt'.format(output_dir, name)).read().split('\n')
        if len(adapter_files[0]) > 0 and not adapter_files[0] == 'None':
            adapter = adapter_files[0].split('/')[-1].split('.fa')[0]
        else:
            adapter_files = list()
            adapter = None

        # For each preprocessing step, a tuple of (prefix, suffix for forward, suffix for reverse)
        prefix2terms = {'[Initial quality assessment]': ('', 'R1', 'R2'),
                        '[Before quality trimming]': (('', 'forward', 'reverse')
                                                      if performed_rrna_removal else ('', adapter + '_forward_paired',
                                                                                      adapter + '_reverse_paired')
                        if adapter is not None else ('', 'R1', 'R2')), '[After quality trimming]': (
                'quality_trimmed_', 'forward_paired', 'reverse_paired')}

        # Initial assessment
        self.info_from_fastqc(output_dir, name, '[Initial quality assessment]', prefix2terms, fastq_columns,
                              initial_file=input_file)

        # After adapter removal
        try:
            if len(adapter_files) > 0:
                self.report.loc[name, '[Adapter removal] adapter files'] = ', '.join(set(adapter_files))
            else:
                self.report.loc[name, '[Adapter removal] adapter files'] = 'None'
        except:
            print('Failed at adapter removal!')
            self.report.to_csv('{}/report.tsv'.format(output_dir), sep='\t')

        # Quality trimming
        self.info_from_fastqc(output_dir, name, '[Before quality trimming]', prefix2terms, fastq_columns,
                              initial_file=input_file)

        self.report.loc[name, '[Quality trimming] Parameters'] = '; '.join([
            file for file in set(open('{}/Preprocess/Trimmomatic/{}_quality_params.txt'.format(output_dir, name)).read(
                ).split('\n')) if len(file) > 2])  # TODO - because '' must be interfering, try to cut the problem at the root before troubles

        self.info_from_fastqc(output_dir, name, '[After quality trimming]', prefix2terms, fastq_columns)

    def set_samples(self, experiments):
        experiments = experiments[['Name', 'Sample']]
        self.report = pd.merge(experiments, self.report, left_on='Name', right_index=True, how='outer')

    def info_from_assembly(self, output_dir, sample):
        print('Retrieving assembly information for sample ' + sample)
        qc_report = pd.read_csv('{}/Assembly/{}/quality_control/report.tsv'.format(
            output_dir, sample), sep='\t', index_col=0).transpose()
        qc_report.index = [sample]

        for col in qc_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Assembly] ' + col] = qc_report.loc[sample][col]

    def info_from_annotation(self, output_dir, sample):
        print('Retrieving annotation information for sample ' + sample)
        sample_report = dict()
        sample_report['# of proteins detected'] = (
            count_on_file('>', '{}/Annotation/{}/fgs.faa'.format(output_dir, sample)))

        sample_report['# of proteins annotated (DIAMOND)'] = (
                len(set(parse_blast('{}/Annotation/{}/aligned.blast'.format(
                    output_dir, sample))['qseqid'])) - count_on_file('*',
                                                                     '{}/Annotation/{}/aligned.blast'.format(output_dir,
                                                                                                             sample)))

        sample_report['# of proteins annotated (reCOGnizer)'] = (
            len(set(parse_blast('{}/Annotation/{}/COG_aligned.blast'.format(
                output_dir, sample))['qseqid'])))

        sample_report = pd.DataFrame.from_dict(sample_report, orient='index').transpose()
        sample_report.index = [sample]

        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Annotation] ' + col] = (
                sample_report.loc[sample][col])

    def info_from_binning(self, output_dir, sample):
        sample_report = dict()
        sample_report['# of bins'] = len(glob.glob(
            '{0}/Binning/{1}/{1}.*.fasta'.format(output_dir, sample)))
        checkm = pd.read_csv('{}/Binning/{}/checkm.tsv'.format(
            output_dir, sample), sep='\t')
        sample_report['# of high-quality drafts'] = ((checkm['Completeness'] > 90)
                                                     & (checkm['Contamination'] < 5)).sum()
        sample_report['# of medium-quality drafts'] = ((checkm['Completeness'] < 90)
                                                       & (checkm['Completeness'] > 50) & (
                                                                   checkm['Contamination'] < 10)).sum()
        sample_report['# of low-quality drafts'] = ((checkm['Completeness'] < 50)
                                                    & (checkm['Contamination'] < 10)).sum()
        sample_report = pd.DataFrame.from_dict(sample_report, orient='index').transpose()
        sample_report.index = [sample]

        for col in sample_report.columns.tolist():
            self.report.loc[self.report['Sample'] == sample, '[Binning] ' + col] = (
                sample_report.loc[sample][col])

    def info_from_alignment(self, output_dir, mt_name):
        self.report.set_index('Name', inplace=True)
        self.report.loc[mt_name, '[Metatranscriptomics] # of reads aligned'] = (
            run_pipe_command("""cat {}/Metatranscriptomics/{}.sam | 
                cut -f 3 | sort | uniq -c | awk \'{{printf(\"%s\\t%s\\n\", $2, $1)}}\' | 
                awk '{{sum+=$2}} END {{print sum}}\'""".format(output_dir, mt_name),
                             output='PIPE').rstrip('\n'))
        self.report.reset_index(inplace=True)

    def info_from_differential_expression(self, output_dir, sample, cutoff=0.01):
        de_results = pd.read_csv('{}/Metatranscriptomics/condition_treated_results.csv'.format(
            output_dir), index_col=0)
        self.report.loc[self.report['Sample'] == sample, '[Gene expression] # of differentially expressed proteins'] = (
            (de_results['padj'] < cutoff).sum())

    def zip_files(self, files, output):
        with ZipFile(output, 'w') as archive:
            for file in files:
                archive.write(file)

    def run(self):
        args = self.get_arguments()

        fastq_columns = open('{}/fastqc_columns.txt'.format(args.lists_directory)).read().split('\n')
        reporter_columns = open('{}/reporter_columns.txt'.format(args.lists_directory)).read().split('\n')

        self.write_technical_report('{}/technical_report.tsv'.format(args.output))

        exps = (
            pd.read_csv(args.experiments, sep='\t') if args.input_format == 'tsv' else pd.read_excel(args.experiments))

        self.initialize_report(reporter_columns)

        for i in exps.index:
            self.info_from_preprocessing(
                args.output, exps.iloc[i]['Name'], exps.iloc[i]['Files'].split(',')[0], fastq_columns,
                performed_rrna_removal=(False if exps.iloc[i]['Data type'] == 'dna' else True))

        self.set_samples(exps)

        self.report.to_csv('report.tsv', sep='\t')
        for sample in set(exps['Sample']):
            self.info_from_assembly(args.output, sample)

            self.info_from_annotation(args.output, sample)

            self.info_from_binning(args.output, sample)

            self.info_from_differential_expression(args.output, sample)

        for mt_name in exps[exps["Data type"] == 'mrna']['Name']:
            self.info_from_alignment(args.output, mt_name)

        self.report.to_excel('{}/MOSCA_General_Report.xlsx'.format(args.output))

        self.zip_files([
            '{}/{}'.format(args.output, filename) for filename in [
                'MOSCA_Protein_Report.xlsx',
                'MOSCA_Entry_Report.xlsx',
                'MOSCA_General_Report.xlsx',
                'technical_report.tsv'
            ]],
            '{}/{}'.format(args.output, 'MOSCA_results.zip'))


if __name__ == '__main__':
    Reporter().run()
