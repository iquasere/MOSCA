# -*- coding: utf-8 -*-
'''
MOSCA Preprocessing package for quality check and removal of
undesired reads by quality, detection of artificial origin
or detection as rRNA

By João Sequeira

March 2017
'''

import shutil
from glob import glob
import os
import pathlib
import time
from subprocess import check_output, DEVNULL
from mosca_tools import run_command, run_pipe_command, parse_fastqc_report, fastqc_name


class Preprocesser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def run_fastqc(self, files, out_dir, extract=True, threads='12'):
        run_command(
            f"fastqc --outdir {out_dir} --threads {threads}{' --extract ' if extract else ' '}{' '.join(files)}")

    def has_adapters(self, fastqc_report):
        data = parse_fastqc_report(fastqc_report)
        if not data['Adapter Content'][0] == 'pass':
            return True
        terms = ['Adapter', 'Illumina', 'Primer']  # Terms that appear in FastQC chapter "overrepresented sequences"
        if not data['Overrepresented sequences'][0] == 'pass':
            i = 0
            while i < len(data['Overrepresented sequences'][1]['Possible Source']):
                for term in terms:
                    if term in data['Overrepresented sequences'][1]['Possible Source'][i]:
                        return True
                i += 1
        return False

    def select_adapters(self, adapters, paired):
        if paired:
            return [adapter for adapter in adapters if 'PE' in adapter]
        return [adapter for adapter in adapters if 'SE' in adapter]

    def download_resources(self, resources_directory):
        if os.path.isfile(f'{resources_directory}/downloaded_timestamp.txt'):
            print(f'{resources_directory}/downloaded_timestamp.txt found! Not downloading resources.')
            return
        print(f'{resources_directory}/downloaded_timestamp.txt not found! Downloading resources.')
        run_pipe_command(
            f"wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz "
            f"-P {resources_directory}/rRNA_databases/")
        run_pipe_command(
            f'tar -xzvf {resources_directory}/rRNA_databases/database.tar.gz -C {resources_directory}/rRNA_databases')
        run_pipe_command(f'cp {self.get_adapters_dir()}/*.fa {resources_directory}/adapters')
        with open(f'{resources_directory}/downloaded_timestamp.txt', 'w') as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))

    def get_adapters_dir(self):
        return glob(f"{'/'.join(check_output('which trimmomatic', shell=True).decode('utf8').split('/')[:-2])}"
                    f"/share/trimmomatic-*/adapters")[0]

    def remove_adapters(self, files, out_dir, name, adapters, threads='12'):
        adapter_contaminated = False
        for file in files:
            folder = fastqc_name(file.split('/')[-1])
            if self.has_adapters(f'{out_dir}/FastQC/{folder}_fastqc/fastqc_data.txt'):
                adapter_contaminated = True
        if not adapter_contaminated:  # No files had adapters
            print('No adapter contamination detected.')
            return 'None'
        for adapter in adapters:  # trim according to each adapter file
            adapter_name = adapter.split('/')[-1].split('.fa')[0]
            output_files = (
                ' '.join([f'{out_dir}/Trimmomatic/noadapters_{name}_{adapter_name}_{fr}_{pu}.fq'
                          for fr in ['forward', 'reverse'] for pu in ['paired', 'unpaired']]) if self.paired else
                f'{out_dir}/Trimmomatic/noadapters_{name}_{adapter_name}.fq')
            run_command(
                f"trimmomatic {'PE' if self.paired else 'SE'} -threads {threads} {' '.join(files)} "
                f"{output_files} ILLUMINACLIP:{adapter}:2:30:10 MINLEN:20")
            self.run_fastqc(
                ([f'{out_dir}/Trimmomatic/noadapters_{name}_{adapter_name}_{fr}_paired.fq'
                    for fr in ['forward', 'reverse']] if self.paired
                    else [f'{out_dir}/Trimmomatic/noadapters_{name}_{adapter_name}.fq']),
                f'{out_dir}/FastQC', threads=threads)
            has_adapters = False
            for file in ([
                f'{out_dir}/FastQC/noadapters_{name}_{adapter_name}_{fr}_paired_fastqc/fastqc_data.txt'
                for fr in ['forward', 'reverse']] if self.paired
                    else [f'{out_dir}/FastQC/noadapters_{name}_{adapter_name}_fastqc/fastqc_data.txt']):
                if self.has_adapters(file):
                    has_adapters = True
            if not has_adapters:
                return adapter  # It's solved, adapters have been removed
        return 'Failed'

    # SortMeRNA - rRNA removal
    def rrna_removal(self, reads, out_dir, name, database, indexes_dir, tmp_dir, threads=12):
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        run_pipe_command(
            f"sortmerna -ref {database} --reads {' --reads '.join(reads)} --idx-dir {indexes_dir} "
            f"--workdir {tmp_dir} --aligned {out_dir}/rrna_{name} --other {out_dir}/norrna_{name} -out2 --fastx "
            f"--paired_in --threads {threads} 1>{out_dir}/{name}_sortmerna.log 2>{out_dir}/{name}_sortmerna.err")
        shutil.rmtree(tmp_dir)
        not_compressed = glob(f'{out_dir}/*.fq')
        if len(not_compressed) > 0:
            for file in not_compressed:
                run_pipe_command(f'gzip {file}')

    def get_crop(self, data):
        data = data['Per base sequence quality'][1]
        i = 0
        while i < len(data):
            if float(data['Lower Quartile'][i]) < 10 or float(data['Median'][i]) < 25:
                crop = data.index[i]
                if '-' in crop:
                    return int(crop.split('-')[0])
                return int(crop)
            else:
                i += 1

    def get_headcrop(self, data):
        def pbsc_value(letter, i):
            return float(data['Per base sequence content'][1][letter][i])
        i = int(len(data['Per base sequence content'][1]) / 2)  # don't want to cut over half of the sequences
        while i > 0 and (abs(pbsc_value('A', i) - pbsc_value('T', i)) < 10 and
                         abs(pbsc_value('G', i) - pbsc_value('C', i)) < 10):
            i -= 1
        headcrop = data['Per base sequence content'][1].index[i]
        if '-' in headcrop:
            return int(headcrop.split('-')[1])
        return int(headcrop) + 1

    # Trimmomatic - removal of low quality regions and short reads
    def quality_trimming(self, files, out_dir, name, threads='12', minlen='100', avgqual='20', type_of_data='dna'):
        fastqc_reports = [
            f"{out_dir}/FastQC/{fastqc_name(file.split('/')[-1])}_fastqc/fastqc_data.txt" for file in files]
        crop = float('inf')
        headcrop = 0
        for report in fastqc_reports:
            data = parse_fastqc_report(report)
            if data['Per base sequence quality'][0] in ['warn', 'fail']:
                parameter = self.get_crop(data)
                if parameter < crop:
                    crop = parameter
            if type_of_data == 'dna':
                if data['Per base sequence content'][0] in ['warn', 'fail']:
                    parameter = self.get_headcrop(data)
                    if parameter > headcrop:
                        headcrop = parameter
        output_files = ' '.join([f'{out_dir}/Trimmomatic/quality_trimmed_{name}_{fr}_{pu}.fq' for fr in [
            'forward', 'reverse'] for pu in ['paired', 'unpaired']]
                                ) if self.paired else f'{out_dir}/Trimmomatic/quality_trimmed_{name}.fq'
        run_command(f"trimmomatic {'PE' if self.paired else 'SE'} -threads {threads} {' '.join(files)} "
                    f"{output_files}{f' CROP:{crop}' if crop < float('inf') else ''}"
                    f"{f' HEADCROP:{headcrop}' if headcrop > 0 else ''} AVGQUAL:{avgqual} MINLEN:{minlen}")

        with open(f'{out_dir}/Trimmomatic/{name}_quality_params.txt', 'a') as f:
            if headcrop > 0:
                f.write(f'HEADCROP:{headcrop}\n')
            if crop < float('inf'):
                f.write(f'CROP:{crop}\n')
            f.write(f'AVGQUAL:{avgqual}\n')
            f.write(f'MINLEN:{minlen}\n')

    # TODO - implement
    def host_sequences_removal(self):
        pass

    def remove_intermediates(self, out_dir):
        for file in (glob(f'{out_dir}/Trimmomatic/noadapters_*.fq') + glob(f'{out_dir}/SortMeRNA/norrna_*.fq.gz')):
            os.remove(file)
            print(f'Removed intermediate file: {file}')

    def run(self):
        reads = snakemake.params.reads.split(',')
        for directory in \
                [f'{snakemake.params.output}/{f}' for f in ['FastQC', 'Trimmomatic', 'SortMeRNA']] + \
                [f'{snakemake.params.resources_directory}{f}' for f in ['/adapters', '/rRNA_databases']]:
            pathlib.Path(directory).mkdir(parents=True, exist_ok=True)
        if snakemake.params.name is None:
            name = fastqc_name(reads[0].split('/')[-1])
            if '_R' in name:
                name = name.split('_R')[0]
        else:
            name = snakemake.params.name
        self.paired = len(reads) > 1

        # First quality check
        self.run_fastqc(reads, f'{snakemake.params.output}/FastQC', threads=snakemake.threads)

        self.download_resources(snakemake.params.resources_directory)

        adapters_dir, rrna_databases_dir = [
            f'{snakemake.params.resources_directory}/{f}' for f in ['adapters', 'rRNA_databases']]

        # Adapter removal
        adapters = self.select_adapters(glob(f'{adapters_dir}/*.fa*'), paired=self.paired)
        print('Available adapter files:\n{}'.format('\n'.join(adapters)))

        adapter_result = self.remove_adapters(reads, snakemake.params.output, name, adapters, threads=snakemake.threads)

        print(f'adapter result:{adapter_result}')
        with open(f'{snakemake.params.output}/Trimmomatic/{name}_adapters.txt', 'w') as f:
            f.write(adapter_result)
        if adapter_result not in ['None', 'Failed']:
            adapter_part = adapter_result.split('/')[-1].split('.fa')[0]
            if self.paired:
                reads = [f'{snakemake.params.output}/Trimmomatic/noadapters_{name}_{adapter_part}_{fr}_paired.fq'
                         for fr in ['forward', 'reverse']]
            else:
                reads = [f'{snakemake.params.output}/Trimmomatic/noadapters_{name}_{adapter_part}.fq']

        # self.host_sequences_removal()

        # rRNA removal
        if snakemake.params.data_type == 'mrna':
            rrna_db = snakemake.params.rrna_db
            rrna_dbs_options = ['fast', 'default', 'sensitive', 'sensitive_with_rfam']
            if rrna_db not in rrna_dbs_options:
                exit(f'rrna_db must be one of {rrna_dbs_options}')
            rrna_database = 'sensitive_db_rfam_seeds' if rrna_db == 'sensitive_with_rfam' else f'{rrna_db}_db'
            rrna_database = f'{rrna_databases_dir}/smr_v4.3_{rrna_database}.fasta'
            self.rrna_removal(
                reads, f'{snakemake.params.output}/SortMeRNA', name, rrna_database, rrna_databases_dir,
                tmp_dir=f'{snakemake.params.output}/SortMeRNA/tmp', threads=snakemake.threads)

            reads = ([f'{snakemake.params.output}/SortMeRNA/norrna_{name}_{fr}.fq.gz' for fr in ['fwd', 'rev']] if
                     self.paired else [f'{snakemake.params.output}/SortMeRNA/norrna_{name}.fq.gz'])

        self.run_fastqc(reads, f'{snakemake.params.output}/FastQC', threads=snakemake.threads)

        self.quality_trimming(
            reads, snakemake.params.output, name, threads=snakemake.threads, avgqual=snakemake.params.avgqual,
            minlen=snakemake.params.mt_minlen if snakemake.params.data_type == 'mrna' else snakemake.params.mg_minlen,
            type_of_data=snakemake.params.data_type)

        reads = ([f'{snakemake.params.output}/Trimmomatic/quality_trimmed_{name}_{fr}_paired.fq'
                  for fr in ['forward', 'reverse']]
                 if self.paired else [f'{snakemake.params.output}/Trimmomatic/quality_trimmed_{name}.fq'])

        self.run_fastqc(reads, f'{snakemake.params.output}/FastQC', threads=snakemake.threads)

        self.remove_intermediates(snakemake.params.output)


if __name__ == '__main__':
    Preprocesser().run()
