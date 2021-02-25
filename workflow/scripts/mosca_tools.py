# -*- coding: utf-8 -*-
'''
General tools for support of MOSCA's functionalities

By João Sequeira

Jun 2017
'''

import glob
import numpy as np
import os
import pandas as pd
import subprocess
import sys
import time


def run_command(bashCommand, output='', mode='w', sep=' ', print_message=True, verbose=True):
    if print_message:
        print('{}{}'.format(bashCommand.replace(sep, ' '), ' > ' + output if output != '' else ''))
    if output == '':
        subprocess.run(bashCommand.split(sep), stdout=sys.stdout if verbose else None,
                       check=True)
    else:
        with open(output, mode) as output_file:
            subprocess.run(bashCommand.split(sep), stdout=output_file)


def run_pipe_command(bashCommand, output='', mode='w', sep=' ', print_message=True):
    if print_message:
        print(bashCommand)
    if output == '':
        subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True,
                                stdout=subprocess.PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            subprocess.Popen(bashCommand, stdin=subprocess.PIPE, shell=True, stdout=output_file).communicate()


def parse_blast(blast):
    result = pd.read_csv(blast, sep='\t', header=None)
    result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                      'gapopen', 'qstart', 'qend', 'sstart', 'send',
                      'evalue', 'bitscore']
    return result


def check_bowtie2_index(index_prefix):
    files = glob.glob(index_prefix + '*.bt2')
    if len(files) < 6:
        return False
    return True


def generate_mg_index(reference, index_prefix):
    run_command('bowtie2-build {} {}'.format(reference, index_prefix))


def align_reads(reads, index_prefix, sam, report, log=None, threads=6):
    run_command('bowtie2 -x {} -1 {} -2 {} -S {} -p {} 1> {} 2> {}'.format(
        index_prefix, reads[0], reads[1], sam, threads, report, log))


def parse_fasta(file):
    print('Parsing', file)
    lines = [line.rstrip('\n') for line in open(file)]
    i = 0
    sequences = dict()
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:]
            sequences[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                sequences[name] += lines[i]
                i += 1
    return sequences


def build_gff_from_contigs(contigs, output, assembler=None):
    gff = pd.DataFrame()
    contigs = parse_fasta(contigs)
    contigs = [[k, v] for k, v in contigs.items()]
    gff["seqid"] = [contig[0] for contig in contigs]
    size = len(contigs)
    gff["source"] = ['.' if assembler is None else assembler] * size
    gff["type"] = ['exon'] * size
    gff["start"] = ['1'] * size
    gff["end"] = [len(contigs[1]) for contigs in contigs]
    gff["score"] = ['.'] * size
    gff["strand"] = ['.'] * size
    gff["phase"] = ['.'] * size
    gff["attributes"] = ['gene_id=' + contig[0] for contig in contigs]
    gff.to_csv(output, sep='\t', index=False, header=False)


def build_gff(blast, output, assembler='metaspades'):
    gff = pd.DataFrame()
    diamond = parse_blast(blast)
    parts = [qid.split('_') for qid in diamond.qseqid]
    preid = [part[1] for part in parts]
    node = 1
    j = 1
    ids = list()
    for i in preid:
        if i == node:
            ids.append('seq' + str(i) + '_' + str(j))
            j += 1
        else:
            node = i
            j = 1
    gff["seqid"] = ['_'.join(part[:-3]) for part in parts]
    size = gff.size
    gff["source"] = ['UniProtKB' for i in range(size)]
    gff["type"] = ['exon' for i in range(size)]
    gff["start"] = [part[-3] for part in parts]
    gff["end"] = [part[-2] for part in parts]
    gff["score"] = diamond.evalue
    gff["strand"] = [part[-1] for part in parts]
    gff["phase"] = ['.' for i in range(size)]
    ids = [ide.split('|')[1] if ide != '*' else ide for ide in diamond.sseqid]
    gff["Name"] = diamond.qseqid
    gff["attributes"] = ['gene_id=' + ids[i] + ';Name=' + diamond.iloc[i]['qseqid'] for i in range(size)]
    del gff["Name"]
    gff.to_csv(output, sep='\t', index=False, header=False)


'''
Input:
    reference: name of contigs file from assembly
    reads: list, [forward reads, reverse reads]
    basename: basename of outputs
    threads: number of threads to use
    blast: name of BLAST file
Output:
    Will generate a bowtie2 index named contigs.replace(.fasta,_index), 
    GFF annotation file named blast.replace(.blast,.gff)
    SAM alignment and READCOUNTS files named output + .sam and .readcounts
'''


def perform_alignment(reference, reads, basename, threads=1, blast=None,
                      blast_unique_ids=True, attribute='gene_id'):
    if not check_bowtie2_index(reference.replace('.fasta', '_index')):
        print('INDEX files not found. Generating new ones')
        generate_mg_index(reference, reference.replace('.fasta', '_index'))
    else:
        print('INDEX was located at ' + reference.replace('.fasta', '_index'))
    if not os.path.isfile(basename + '.log'):
        align_reads(reads, reference.replace('.fasta', '_index'), basename + '.sam',
                    basename + '_bowtie2_report.txt', log=basename + '.log', threads=threads)
    else:
        print('{}.log was found!'.format(basename))

    if blast is None:
        if not os.path.isfile(reference.replace('.fasta', '.gff')):
            print('GFF file not found at {}. Generating a new one.'.format(reference.replace('.fasta', '.gff')))
            build_gff_from_contigs(reference, reference.replace('.fasta', '.gff'))
        else:
            print('GFF file was located at {}'.format(reference.replace('.fasta', '.gff')))
    else:
        if not os.path.isfile(blast.replace('.blast', '.gff')):
            print('GFF file not found at {}. Generating a new one.'.format(reference.replace('.blast', '.gff')))
            if blast_unique_ids:
                build_gff(blast, blast.replace('.blast', '.gff'))
        else:
            print('GFF file was located at ' + blast.replace('.blast', '.gff'))

    run_command('htseq-count -i {0} -c {1}.readcounts -n {2} {1}.sam {3}{4}'.format(
        attribute, basename, threads, (reference.replace('.fasta', '.gff') if blast is None else
                                       blast.replace('.blast', '.gff')),
        ('' if blast is not None else ' --stranded=no')))


def fastq2fasta(fastq, output):
    run_command("paste - - - - < {} | cut -f 1,2 | sed 's/^@/>/' | tr \"\t" "\n\" > {}".format(
        fastq, output))


def timed_message(message):
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ': ' + message)


def normalize_mg_readcounts_by_size(readcounts, contigs):
    run_pipe_command('head -n -5 {}'.format(readcounts),
                     # I still don't know how to pipe two things together for the join command, when I do I'll be able to merge this two commands
                     output=readcounts.replace('.readcounts', '_no_tail.readcounts'))
    run_pipe_command(
        "seqkit fx2tab {} | sort | awk '{{print $1\"\\t\"length($2)}}' | join - {} | awk '{{print $1\"\\t\"$3/$2}}'".format(
            contigs, readcounts.replace('.readcounts', '_no_tail.readcounts')),
        output=readcounts.replace('.readcounts', '_normalized.readcounts'))


'''
Input: 
    data: pd.DataFrame - Protein Report on the making
    readcounts: readcounts file from quantification with htseq-count
    origin_of_data: 'metagenomics', 'metatranscriptomics' defines column of joining
    name: name of column to add to relation
    readcounts_has_tail: boolean, if file is htseq-count expression matrix with
    those last 5 lines of general information
Output: 
    'data' df will receive additional column with abundance/expression 
    information. This column will be named 'name'
'''


def add_abundance(data, readcounts, name, origin_of_data='metagenomics', readcounts_has_tail=True):
    readcounts = pd.read_csv(readcounts, sep='\t', header=None, names=['qseqid', name])
    readcounts[name] = readcounts[name].fillna(value=0)
    if readcounts_has_tail:
        readcounts = readcounts[:-5]  # remove those last resume lines from htseq-count
    if origin_of_data == 'metagenomics':
        readcounts['Contig'] = [qseqid.split('_')[1] for qseqid in readcounts['qseqid']]
        del readcounts['qseqid']
        return pd.merge(data, readcounts, on='Contig', how='left')
    elif origin_of_data == 'metatranscriptomics':
        return pd.merge(data, readcounts, on='qseqid', how='left')
    elif origin_of_data == 'compomics':
        pass


def multi_sheet_excel(output, data, sheet_name='Sheet', lines=1000000, index=False):
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    if len(data) < lines:
        data.to_excel(writer, sheet_name='{}'.format(sheet_name), index=index)
    else:
        for i in range(0, len(data), lines):
            j = min(i + lines, len(data))
            data.iloc[i:(i + lines)].to_excel(writer, sheet_name='{} ({})'.format(sheet_name, j), index=index)
    writer.save()


'''
Input:
    joined: name of final TSV file outputed by MOSCA
    columns: name of columns to normalize (don't put them all at once, one
    normalization at a time!)
    output: name of file to output
Output:
    A TXT file will be generated with a single column of numbers, each one
    corresponding to the normalization factor of each sample by order of
    column in the readcounts file
'''


def normalize_readcounts(joined, columns, method='TMM', rscript_folder=''):
    working_dir = '/'.join(joined.split('/')[:-1])
    info = pd.read_csv(joined, sep='\t')
    info[columns] = info[columns].fillna(value=0)
    info[columns].to_csv(working_dir + '/to_normalize.tsv', sep='\t', index=False)
    print('Normalizing {} on columns {}'.format(joined, ','.join(columns)))
    run_command(
        '{3}Rscript {1}/normalization.R --readcounts {0}/to_normalize.tsv --output {0}/normalization_factors.txt -m {2}'.format(
            working_dir, sys.path[0], method, rscript_folder))
    factors = open(working_dir + '/normalization_factors.txt').read().split('\n')[:-1]  # \n always last element

    for i in range(len(columns)):
        info[columns[i] + '_normalized'] = info[columns[i]] * float(factors[i])
    return info


def timed_message(message=None):
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ': ' + message)


def expand_by_list_column(self, df, column='Pathway'):
    lens = [len(item) for item in df[column]]
    dictionary = dict()
    for column in df.columns:
        dictionary[column] = np.repeat(df[column].values, lens)
    dictionary[column] = np.concatenate(df[column].values)
    return pd.DataFrame(dictionary)


def parse_fastqc_report(filename):
    data = dict()
    file = open(filename).read().split('\n')
    i = 1
    while i < len(file):
        if file[i].startswith('>>') and file[i] != '>>END_MODULE':
            name, flag = file[i][2:].split('\t')[0], file[i][2:].split('\t')[1]
            if name == 'Sequence Duplication Levels':
                i += 1
            i += 1
            if file[i] == '>>END_MODULE':
                data[name] = (flag, pd.DataFrame())
            else:
                labels = file[i][1:].split('\t')
                i += 1
                partial_data = np.array(labels)
                while i < len(file) and not file[i].startswith('>>'):
                    partial_data = np.append(partial_data, file[i].split('\t'))
                    i += 1
                partial_data = np.reshape(partial_data, (int(partial_data.size / len(labels)), len(labels)))
                data[name] = (flag, pd.DataFrame(data=partial_data[1:, 1:], index=partial_data[1:, 0],
                                                 columns=partial_data[0, 1:]))
        i += 1
    return data


'''
Input:
    bashCommand: str - the command to retrieve the output from
    shell: bool - True if using some shell tool like grep or awk
Output:
    Number of occurrences of character on file
'''


def count_on_file(expression, file, compressed=False):
    return int(subprocess.check_output("{} -c '{}' {}".format(
        'zgrep' if compressed else 'grep', expression, file), shell=True))


def sort_alphanumeric(alphanumeric_list):
    return sorted(alphanumeric_list, key=lambda item: (int(item.partition(' ')[0])
                                                       if item[0].isdigit() else float('inf'), item))
