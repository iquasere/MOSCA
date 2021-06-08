# -*- coding: utf-8 -*-
'''
MOSCA Preprocessing package for quality check and removal of
undesired reads by quality, detection of artificial origin
or detection as rRNA

By João Sequeira

March 2017
'''

import argparse
import glob
import multiprocessing
import os
import sys
import pathlib
import time
from mosca_tools import run_command, run_pipe_command, parse_fastqc_report


class Preprocesser:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="MOSCA preprocessing")
        parser.add_argument("-i", "--input", type=str, required=True,
                            help="File(s) for preprocessing")
        parser.add_argument("-t", "--threads", type=str,
                            default=str(multiprocessing.cpu_count() - 2),
                            help="Number of threads to use. Default is number of CPUs available minus 2.")
        parser.add_argument("-d", "--data", type=str, choices=["dna", "mrna"])
        parser.add_argument("-o", "--output", type=str, help="Output directory"),
        parser.add_argument("-adaptdir", "--adapters-directory", type=str,
                            help="Directory with adapter files for Trimmomatic",
                            default=os.path.expanduser('~/illumina_adapters'))
        parser.add_argument("-rrnadbs", "--rrna-databases-directory", type=str,
                            help="Directory with rRNA databases for SortMeRNA",
                            default=os.path.expanduser('~/rRNA_databases'))
        parser.add_argument("-rd", "--resources-directory", type=str,
                            help="Directory with resources for SortMeRNA and Trimmomatic")
        parser.add_argument("-n", "--name", type=str,
                            help="Name attributed to the data inputted. Will be added to basename of output files")
        parser.add_argument("--minlen", default=100, help='Minimum length of reads to keep')
        parser.add_argument("--avgqual", default=20, help='Average quality of reads to keep')

        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        args.input = args.input.split(',')
        return args

    def run_fastqc(self, files, out_dir, extract=True, threads='12'):
        run_command('fastqc --outdir {} --threads {}{}{}'.format(
            out_dir, threads, ' --extract ' if extract else ' ', ' '.join(files)))

    def has_adapters(self, fastqc_report):
        data = parse_fastqc_report(fastqc_report)
        if not data['Adapter Content'][0] == 'pass':
            return True
        terms_list = [
            'Adapter', 'Illumina', 'Primer']  # Terms that appear in FastQC chapter "overrepresented sequences"
        if not data['Overrepresented sequences'][0] == 'pass':
            i = 0
            while i < len(data['Overrepresented sequences'][1]['Possible Source']):
                for term in terms_list:
                    if term in data['Overrepresented sequences'][1]['Possible Source'][i]:
                        return True
                i += 1
        return False

    def select_adapters(self, adapters, paired):
        if paired:
            return [adapter for adapter in adapters if 'PE' in adapter]
        else:
            return [adapter for adapter in adapters if 'SE' in adapter]

    def fastqc_name(self, filename):
        return filename.replace("stdin:", "").replace(".gz", "").replace(".bz2", "").replace(".txt", "").replace(
            ".fastq", "").replace(".fq", "").replace(".csfastq", "").replace(".sam", "").replace(".bam", "")

    def remove_fa_end(self, filename):
        if 'fasta' in filename:
            return filename.split('.fasta')[0]
        else:
            return filename.split('.fa')[0]

    def download_resources(self, resources_directory):
        if not os.path.isfile(f'{resources_directory}/downloaded_timestamp.txt'):
            run_command(f'bash {sys.path[0]}/download_resources.sh {resources_directory}')
            with open(f'{resources_directory}/downloaded_timestamp.txt', 'w') as f:
                f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        else:
            print(f'{resources_directory}/downloaded_timestamp.txt found! Not downloading resources.')

    '''
    Input:
    Output:
        If files had adapters, and those were successfully removed, returns adapter 
        filename used. If adapters could not be removed, returns 'Failed'. If no
        adapter presence was detected from the start, returns 'None'
    '''

    def remove_adapters(self, files, out_dir, name, adapters, threads='12'):
        adapter_contaminated = False
        for file in files:
            folder = self.fastqc_name(file.split('/')[-1])
            if self.has_adapters(f'{out_dir}/FastQC/{folder}_fastqc/fastqc_data.txt'):
                adapter_contaminated = True
        if not adapter_contaminated:  # No files had adapters
            print('No adapter contamination detected.')
            return 'None'
        for adapter in adapters:  # trim according to each adapter file
            adapter_name = self.remove_fa_end(adapter.split('/')[-1])

            run_command('trimmomatic {} -threads {} {} {} ILLUMINACLIP:{}:2:30:10 MINLEN:20'.format(
                'PE' if self.paired else 'SE', threads, ' '.join(files),
                ' '.join(['{}/Trimmomatic/after_adapter_removal_{}_{}_{}_{}.fq'.format(
                    out_dir, name, adapter_name, fr, pu) for fr in ['forward', 'reverse']
                    for pu in ['paired', 'unpaired']]) if self.paired else
                '{}/Trimmomatic/after_adapter_removal_{}_{}.fq'.format(
                    out_dir, name, adapter_name), adapter))

            self.run_fastqc(
                ([f'{out_dir}/Trimmomatic/after_adapter_removal_{name}_{adapter_name}_{fr}_paired.fq'
                  for fr in ['forward', 'reverse']] if self.paired
                 else [f'{out_dir}/Trimmomatic/after_adapter_removal_{name}_{adapter_name}.fq']),
                f'{out_dir}/FastQC', threads=threads)

            has_adapters = False
            for file in ([
                f'{out_dir}/FastQC/after_adapter_removal_{name}_{adapter_name}_{fr}_paired_fastqc/fastqc_data.txt'
                for fr in ['forward', 'reverse']] if self.paired
            else [f'{out_dir}/FastQC/after_adapter_removal_{name}_{adapter_name}_fastqc/fastqc_data.txt']):
                if self.has_adapters(file):
                    has_adapters = True
            if not has_adapters:
                return adapter  # It's solved, adapters have been removed
        return 'Failed'

    def index_rrna_database(self, database):
        for termination in ['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats']:
            if not os.path.isfile(database.replace('.fasta', '.idx.{}'.format(termination))):
                run_command('indexdb_rna --ref {},{}'.format(database, database.replace('fasta', 'idx')))

    def merge_pe(self, forward, reverse, interleaved):
        run_command(f'bash {sys.path[0]}/merge-paired-reads.sh {forward} {reverse} {interleaved}')

    def unmerge_pe(self, interleaved, forward, reverse):
        run_command(f'bash {sys.path[0]}/unmerge-paired-reads.sh {interleaved} {forward} {reverse}')

    '''
    Input:
        filename: str - filename of FastQ file with irregular reads (less than
                        4 lines per read)
    Output:
        Irregular reads will be removed
    '''

    def remove_messed_reads(self, filename):
        run_pipe_command("awk 'BEGIN {{RS=\"@\"; FS=\"\\n\"}}{{if (NF == 5) print \"@\" substr($0,1,length-1)}}' {0} "
                         "> {0}.temp".format(filename))
        os.rename(f"{filename}.temp", filename)

    # correct number of reads per file - if unequal number of reads from forward to reverse file, it will be corrected by separation name/1,2
    # from www.biostars.org/p/6925/#6928
    def remove_orphans(self, forward, reverse, out_dir):
        run_pipe_command(
            """awk '{{printf $0;getline; printf \"\\t\"$0;getline;getline;print \"\\t\"$0}}' {} | sort -T.""".format(
                forward), output="{}/read1.txt".format(out_dir))

        run_pipe_command(
            """awk '{{printf $0;getline; printf \"\\t\"$0;getline;getline;print \"\\t\"$0}}' {} | sort -T.""".format(
                forward), output="{}/read2.txt".format(out_dir))

        run_pipe_command(
            """join {} | awk '{{print $1\" \"$2\"\\n\"$3\"\\n+\\n\"$4 > \"{}\";print $1\" \"$5\"\\n\"$6\"\\n+\\n\"$7 > \"{}\"}}'""".format(
                ' '.join(["{}/{}".format(out_dir, fr) for fr in ['read1.txt', 'read2.txt']]), forward, reverse))

    # SortMeRNA - rRNA removal
    def rrna_removal(self, files, out_dir, name, databases, threads='12'):
        for database in databases:
            self.index_rrna_database(database)

        for i in range(len(files)):
            if '.gz' in files[i]:
                run_command(f'gunzip {files[i]}')
                files[i] = files[i].rstrip('.gz')

        if self.paired:
            tool_input = f'{out_dir}/{name}_interleaved.fastq'
            self.merge_pe(files[0], files[1], tool_input)
        else:
            tool_input = files[0]

        run_command(
            'sortmerna --ref {0} --reads {1} --aligned {2}/after_rrna_removal_{3}_accepted --fastx '
            '--other {2}/after_rrna_removal_{3}_rejected -a {4}{5}'.format(
                ':'.join(['{0}.fasta,{0}.idx'.format(self.remove_fa_end(database)) for database in databases]),
                tool_input, out_dir, name, threads, ' --paired_out' if len(files) > 1 else ''))

        if self.paired:
            self.unmerge_pe(f'{out_dir}/after_rrna_removal_{name}_rejected.fastq',
                            f'{out_dir}/after_rrna_removal_{name}_forward.fastq',
                            f'{out_dir}/after_rrna_removal_{name}_reverse.fastq')

            for fr in ['forward', 'reverse']:
                self.remove_messed_reads(f'{out_dir}/after_rrna_removal_{name}_{fr}.fastq')

            self.remove_orphans(f'{out_dir}/after_rrna_removal_{name}_forward.fastq',
                                f'{out_dir}/after_rrna_removal_{name}_reverse.fastq',
                                out_dir)

    def get_crop(self, data):
        i = 0
        while i < len(data['Per base sequence quality'][1]):
            if (float(data['Per base sequence quality'][1]['Lower Quartile'][i]) < 10
                    or float(data['Per base sequence quality'][1]['Median'][i]) < 25):
                crop = data['Per base sequence quality'][1].index[i]
                if '-' in crop:
                    return int(crop.split('-')[0])
                return int(crop)
            else:
                i += 1

    def get_headcrop(self, data):
        i = int(len(data['Per base sequence content'][1]) / 2)  # don't want to cut over half of the sequences
        while i > 0 and (abs(float(data['Per base sequence content'][1]['A'][i]) -
                             float(data['Per base sequence content'][1]['T'][i])) < 10
                         and abs(float(data['Per base sequence content'][1]['G'][i]) -
                                 float(data['Per base sequence content'][1]['C'][i])) < 10):
            i -= 1
        headcrop = data['Per base sequence content'][1].index[i]
        if '-' in headcrop:
            return int(headcrop.split('-')[1])
        return int(headcrop) + 1

    # Trimmomatic - removal of low quality regions and short reads
    def quality_trimming(self, files, out_dir, name, threads='12', minlen='100', avgqual='20', type_of_data='dna'):
        fastqc_reports = [f"{out_dir}/FastQC/{self.fastqc_name(file.split('/')[-1])}_fastqc/fastqc_data.txt"
                          for file in files]

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

        files = ' '.join([f'{out_dir}/Trimmomatic/quality_trimmed_{name}_{fr}_{pu}.fq' for fr in ['forward', 'reverse']
                         for pu in ['paired', 'unpaired']]) if self.paired else \
            f'{out_dir}/Trimmomatic/quality_trimmed_{name}.fq'
        run_command(f"trimmomatic {'PE' if self.paired else 'SE'} -threads {threads} {' '.join(files)}"
                    f"{files}{f' CROP:{crop}' if crop < float('inf') else ''}"
                    f"{f' HEADCROP:{headcrop}' if headcrop > 0 else ''}"
                    f" AVGQUAL:{avgqual} MINLEN:{minlen}")

        with open(f'{out_dir}/Trimmomatic/{name}_quality_params.txt', 'a') as f:
            if headcrop > 0:
                f.write(f'HEADCROP:{headcrop}\n')
            if crop < float('inf'):
                f.write(f'CROP:{crop}\n')
            f.write(f'AVGQUAL{avgqual}\n')
            f.write(f'MINLEN:{minlen}\n')

    # TODO - implement
    def host_sequences_removal(self):
        pass

    def remove_intermediates(self, out_dir):
        for file in (glob.glob(f'{out_dir}/Trimmomatic/after_adapter_removal_*.fq') +
                     glob.glob(f'{out_dir}/SortMeRNA/after_rrna_removal_*.fastq') +
                     glob.glob(f'{out_dir}/SortMeRNA/read*.txt')):
            os.remove(file)
            print(f'Removed intermediate file: {file}')

    def run(self):
        args = self.get_arguments()
        original_input = args.input

        for directory in ['FastQC', 'Trimmomatic', 'SortMeRNA']:
            pathlib.Path(f'{args.output}/{directory}').mkdir(parents=True, exist_ok=True)

        name = args.name
        if args.name is None:
            name = self.fastqc_name(args.input[0].split('/')[-1])
            if '_R' in name:
                name = name.split('_R')[0]

        self.paired = True if len(args.input) > 1 else False

        # First quality check
        self.run_fastqc(args.input, f'{args.output}/FastQC', threads=args.threads)

        if args.resources_directory is not None:
            self.download_resources(args.resources_directory)
            args.adapters_directory = f'{args.resources_directory}/adapters'
            args.rrna_databases_directory = f'{args.resources_directory}/rRNA_databases'

        # Adapter removal
        adapters = self.select_adapters(glob.glob(f'{args.adapters_directory}/*.fa*'), paired=self.paired)
        print('Available adapter files:\n{}'.format('\n'.join(adapters)))

        adapter_result = self.remove_adapters(args.input, args.output, name, adapters, threads=args.threads)

        print('adapter result:{}'.format(adapter_result))
        with open(f'{args.output}/Trimmomatic/{name}_adapters.txt', 'w') as f:
            f.write(adapter_result)

        if adapter_result not in ['None', 'Failed']:
            adapter_part = self.remove_fa_end(adapter_result.split('/')[-1])
            if self.paired:
                args.input = [f'{args.output}/Trimmomatic/after_adapter_removal_{name}_{adapter_part}_{fr}_paired.fq'
                              for fr in ['forward', 'reverse']]
            else:
                args.input = [f'{args.output}/Trimmomatic/after_adapter_removal_{name}_{adapter_part}.fq']

        # self.host_sequences_removal()

        # rRNA removal
        if args.data == 'mrna':
            rrna_databases = glob.glob('{}/*.fa*'.format(args.rrna_databases_directory))
            self.rrna_removal(args.input, f'{args.output}/SortMeRNA', name, rrna_databases, threads=args.threads)

            if self.paired:
                args.input = [
                    f'{args.output}/SortMeRNA/after_rrna_removal_{name}_{fr}.fastq' for fr in ['forward', 'reverse']]
            else:
                args.input = [f'{args.output}/SortMeRNA/after_rrna_removal_{name}_rejected.fastq']

        self.run_fastqc(args.input, f'{args.output}/FastQC', threads=args.threads)

        self.quality_trimming(args.input, args.output, name, threads=args.threads, avgqual=args.avgqual,
                              minlen=args.minlen, type_of_data=args.data)

        if self.paired:
            args.input = [
                f'{args.output}/Trimmomatic/quality_trimmed_{name}_{fr}_paired.fq' for fr in ['forward', 'reverse']]
        else:
            args.input = [f'{args.output}/Trimmomatic/quality_trimmed_{name}.fq']

        self.run_fastqc(args.input, f'{args.output}/FastQC', threads=args.threads)

        self.remove_intermediates(args.output)


if __name__ == '__main__':
    Preprocesser().run()
