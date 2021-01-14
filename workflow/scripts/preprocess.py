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
from mosca_tools import run_command, run_pipe_command, parse_fastqc


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

        args = parser.parse_args()

        args.output = args.output.rstrip('/')
        args.input = args.input.split(',')
        return args

    def run_fastqc(self, files, out_dir, extract=True, threads='12'):
        run_command('fastqc --outdir {} --threads {}{}{}'.format(out_dir, threads,
                                                                 ' --extract ' if extract else ' ', ' '.join(files)))

    def has_adapters(self, fastqc_report):
        data = parse_fastqc(fastqc_report)
        terms_list = ['Adapter', 'Illumina',
                      'Primer']  # Terms that appear in FastQC chapter "overrepresented sequences"
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

    def remove_fq_end(self, filename):
        if 'fastq' in filename:
            return filename.split('.fastq')[0]
        else:
            return filename.split('.fq')[0]

    def remove_fa_end(self, filename):
        if 'fasta' in filename:
            return filename.split('.fasta')[0]
        else:
            return filename.split('.fa')[0]

    def download_resources(self, resources_directory):
        if not os.path.isfile('{}/downloaded_timestamp.txt'.format(resources_directory)):
            run_command('bash {}/download_resources.sh {}'.format(sys.path[0], resources_directory))
            with open('{}/downloaded_timestamp.txt'.format(resources_directory), 'w') as f:
                f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        else:
            print('{}/downloaded_timestamp.txt found! Not downloading resources.'.format(resources_directory))

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
            folder = self.remove_fq_end(file.split('/')[-1])  # This is how FastQC produces folder names
            if self.has_adapters('{}/FastQC/{}_fastqc/fastqc_data.txt'.format(
                    out_dir, folder)):
                adapter_contaminated = True
        if not adapter_contaminated:  # No files had adapters
            print('No adapter contamination detected.')
            return 'None'
        for adapter in adapters:  # trim according to each adapter file
            adapter_name = self.remove_fa_end(adapter.split('/')[-1])
            run_command('trimmomatic {} -threads {} {} {} ILLUMINACLIP:{}:2:30:10'.format(
                'PE' if self.paired else 'SE', threads, ' '.join(files),
                ' '.join(['{}/Trimmomatic/after_adapter_removal_{}_{}_{}_{}.fq'.format(
                    out_dir, name, adapter_name, fr, pu) for fr in ['forward', 'reverse']
                    for pu in ['paired', 'unpaired']]) if self.paired else
                '{}/Trimmomatic/after_adapter_removal_{}_{}.fq'.format(
                    out_dir, name, adapter_name), adapter))

            self.run_fastqc((['{}/Trimmomatic/after_adapter_removal_{}_{}_{}_paired.fq'.format(
                out_dir, name, adapter_name, fr) for fr in ['forward', 'reverse']] if self.paired
                             else ['{}/Trimmomatic/after_adapter_removal_{}_{}.fq'.format(
                out_dir, name, adapter_name)]),
                            '{}/FastQC'.format(out_dir), threads=threads)

            has_adapters = False
            for file in ['{}/FastQC/{}_{}_{}_paired_fastqc/fastqc_data.txt'.format(
                    out_dir, name, adapter_name, fr) for fr in ['forward', 'reverse']]:
                if has_adapters(file):
                    has_adapters = True
            if not has_adapters:
                return adapter  # It's solved, adapters have been removed
        return 'Failed'

    def index_rrna_database(self, database):
        for termination in ['bursttrie_0.dat', 'kmer_0.dat', 'pos_0.dat', 'stats']:
            if not os.path.isfile(database.replace('.fasta', '.idx.{}'.format(termination))):
                run_command('indexdb_rna --ref {},{}'.format(database, database.replace('fasta', 'idx')))

    def merge_pe(self, forward, reverse, interleaved):
        run_command('bash {}/merge-paired-reads.sh {} {} {}'.format(
            sys.path[0], forward, reverse, interleaved))

    def unmerge_pe(self, interleaved, forward, reverse):
        run_command('bash {}/unmerge-paired-reads.sh {} {} {}'.format(
            sys.path[0], interleaved, forward, reverse))

    '''
    Input:
        filename: str - filename of FastQ file with irregular reads (less than
                        4 lines per read)
    Output:
        Irregular reads will be removed
    '''

    def remove_messed_reads(self, filename):
        run_pipe_command("""awk 'BEGIN {{RS=\"@\"; FS=\"\\n\"}}{{if (NF == 5)
        "print \"@\" substr($0,1,length-1)}}' {0} > {0}.temp""".format(filename))
        os.rename("{}.temp".format(filename), filename)

    # correct number of reads per file - if unequal number of reads from forward to reverse file, it will be corrected by separation name/1,2
    # from www.biostars.org/p/6925/#6928
    def remove_orphans(self, forward, reverse, out_dir):
        run_pipe_command("""awk '{{printf substr($0,1,length-2);getline;
            printf \"\\t\"$0;getline;getline;print \"\\t\"$0}}' {} | sort -T. > 
            {}/SortMeRNA/read1.txt""".format(forward, out_dir))

        run_pipe_command("""awk '{{printf substr($0,1,length-2);getline;
            printf \"\\t\"$0;getline;getline;print \"\\t\"$0}}' {} | sort -T. > 
            {}/SortMeRNA/read2.txt""".format(forward, out_dir))

        run_pipe_command("""join {} | awk '{{print $1\" \"$2\"\\n\"$3\"\\n+\\n\"$4 > 
            \"{}\";print $1\" \"$5\"\\n\"$6\"\\n+\\n\"$7 > \"{}\"}}'""".format(
            ' '.join(["{}/SortMeRNA/{}".format(out_dir, fr)
                      for fr in ['read1.txt', 'read2.txt']]), forward, reverse))

        for file in ["{}/SortMeRNA/read{}.txt".format(out_dir, number)
                     for number in ['1', '2']]:
            os.remove(file)

    # SortMeRNA - rRNA removal
    def rrna_removal(self, files, out_dir, name, databases, threads='12'):
        for database in databases:
            self.index_rrna_database(database)

        for i in range(len(files)):
            if '.gz' in files[i]:
                run_command('gunzip {}'.format(files[i]))
                files[i] = files[i].rstrip('.gz')

        if self.paired:
            tool_input = '{}/{}_interleaved.fastq'.format(out_dir, name)
            self.merge_pe(files[0], files[1], tool_input)

        else:
            tool_input = files[0]

        run_command(
            'sortmerna --ref {0} --reads {1} --aligned {2}/{3}_accepted --fastx --other {2}/{3}_rejected -a {4}{5}'.format(
                ':'.join(['{0}.fasta,{0}.idx'.format(self.remove_fa_end(database)) for database in databases]),
                tool_input, out_dir, name, threads, ' --paired_out' if len(files) > 1 else ''))

        if self.paired:
            self.unmerge_pe(tool_input, '{}/{}_forward.fastq'.format(out_dir, name),
                            '{}/{}_reverse.fastq'.format(out_dir, name))

        # TODO - check if this is still needed, using awk
        '''
        for fr in ['forward', 'reverse']:
                self.remove_messed_reads('{}_{}.fastq'.format(basename, fr))
            
            self.remove_orphans(basename + '_forward.fastq', 
                               basename + '_reverse.fastq')
        '''

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
    def quality_trimming(self, files, out_dir, name, threads='12'):
        fastqc_reports = ['{}/FastQC/{}_fastqc/fastqc_data.txt'.format(
            out_dir, self.remove_fq_end(file.split('/')[-1])) for file in files]

        crop = float('inf')
        headcrop = 0
        for report in fastqc_reports:
            data = parse_fastqc(report)
            if data['Per base sequence quality'][0] in ['warn', 'fail']:
                parameter = self.get_crop(data)
                if parameter < crop:
                    crop = parameter
            if data['Per base sequence content'][0] in ['warn', 'fail']:
                parameter = self.get_headcrop(data)
                if parameter > headcrop:
                    headcrop = parameter

        run_command('trimmomatic {} -threads {} {} {}{}{} AVGQUAL:20 MINLEN:100'.format(
            'PE' if self.paired else 'SE', threads, ' '.join(files),
            ' '.join('{}/Trimmomatic/quality_trimmed_{}_{}_{}.fq'.format(
                out_dir, name, fr, pu) for fr in ['forward', 'reverse'] for pu in ['paired', 'unpaired'])
            if self.paired else '{}/Trimmomatic/quality_trimmed_{}.fq'.format(
                out_dir, name),
            ' CROP:{}'.format(crop) if crop < float('inf') else '',
            ' HEADCROP:{}'.format(headcrop) if headcrop > 0 else ''))

        with open('{}/Trimmomatic/{}_quality_params.txt'.format(
                out_dir, name), 'a') as f:
            if headcrop > 0: f.write('HEADCROP:{}\n'.format(headcrop))
            if crop < float('inf'): f.write('CROP:{}\n'.format(crop))
            f.write('AVGQUAL{}\n'.format(20))
            f.write('MINLEN:{}\n'.format(100))

    # TODO - implement
    def host_sequences_removal(self):
        pass

    def run(self):
        args = self.get_arguments()

        for directory in ['FastQC', 'Trimmomatic', 'SortMeRNA']:
            pathlib.Path('{}/{}'.format(args.output, directory)).mkdir(parents=True, exist_ok=True)

        name = args.name
        if args.name is None:
            name = self.remove_fq_end(args.input[0].split('/')[-1])
            if '_R' in name: name = name.split('_R')[0]

        self.paired = True if len(args.input) > 1 else False

        # First quality check
        self.run_fastqc(args.input, '{}/FastQC'.format(args.output), threads=args.threads)

        if args.resources_directory is not None:
            self.download_resources(args.resources_directory)
            args.adapters_directory = '{}/adapters'.format(args.resources_directory)
            args.rrna_databases_directory = '{}/rRNA_databases'.format(args.resources_directory)

        # Adapter removal
        adapters = self.select_adapters(glob.glob('{}/*.fa*'.format(args.adapters_directory)), paired=self.paired)
        print('Available adapter files:\n{}'.format('\n'.join(adapters)))

        adapter_result = self.remove_adapters(args.input, args.output, name, adapters, threads=args.threads)
        with open('{}/Trimmomatic/{}_adapters.txt'.format(
                args.output, name), 'w') as f:
            f.write(adapter_result)

        if adapter_result not in ['None', 'Failed']:
            adapter_part = self.remove_fa_end(adapter_result[0].split('/')[-1])
            if self.paired:
                args.input = ['{}/Trimmomatic/after_adapter_removal_{}_{}_{}_paired.fq'.format(
                    args.output, name, adapter_part, fr) for fr in ['forward', 'reverse']]
            else:
                args.input = ['{}/Trimmomatic/after_adapter_removal_{}_{}.fq'.format(
                    args.output, name, adapter_part)]

        # self.host_sequences_removal()

        # rRNA removal
        if args.data == 'mrna':
            rrna_databases = glob.glob('{}/*.fa*'.format(args.rrna_databases_directory))
            self.rrna_removal(args.input, args.output + '/SortMeRNA', name, rrna_databases, threads=args.threads)

            if self.paired:
                args.input = ['{}/SortMeRNA/{}_{}.fastq'.format(
                    args.output, name, fr) for fr in ['forward', 'reverse']]
            else:
                args.input = ['{}/SortMeRNA/{}_rejected.fastq'.format(
                    args.output, name)]

        self.run_fastqc(args.input, '{}/FastQC'.format(args.output), threads=args.threads)

        self.quality_trimming(args.input, args.output, name, threads=args.threads)

        if self.paired:
            args.input = ['{}/Trimmomatic/quality_trimmed_{}_{}_paired.fq'.format(args.output, name,
                                                                                  fr) for fr in ['forward', 'reverse']]
        else:
            args.input = ['{}/Trimmomatic/quality_trimmed_{}.fq'.format(args.output, name)]

        self.run_fastqc(args.input, '{}/FastQC'.format(args.output), threads=args.threads)


if __name__ == '__main__':
    Preprocesser().run()
