# -*- coding: utf-8 -*-
'''
General tools for support of MOSCA's functionalities

By João Sequeira

Jun 2017
'''

from tqdm import tqdm
import pandas as pd, subprocess, glob, os, gzip, time, numpy as np, shutil, sys

class MoscaTools:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

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
    def add_abundance(self, data, readcounts, name, origin_of_data = 'metagenomics', 
                      readcounts_has_tail = True):
        readcounts = pd.read_csv(readcounts, sep = '\t', header = None, 
                                 names = ['qseqid', name])
        readcounts[name] = readcounts[name].fillna(value = 0)
        if readcounts_has_tail: readcounts = readcounts[:-5]                    # remove those last resume lines from htseq-count
        if origin_of_data == 'metagenomics':
            readcounts['Contig'] = [qseqid.split('_')[1] for qseqid in readcounts['qseqid']]
            del readcounts['qseqid']
            return pd.merge(data, readcounts, on = 'Contig', how = 'left')
        elif origin_of_data == 'metatranscriptomics':
            return pd.merge(data, readcounts, on = 'qseqid', how = 'left')
        elif origin_of_data == 'compomics':
            pass
    
    '''
    Input:
        fastq: FASTQ reads filename to convert to FASTA
        output: name of FASTA file to produce
    Output:
        A FASTA version of the input will be created named 'output'
    '''
    def fastq2fasta(self, fastq, output):
        self.run_command("paste - - - - < " + fastq + "| cut -f 1,2 | sed " + 
                         "'s/^@/>/' | tr \"\t" "\n\" > " + output)
    
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
    def normalize_readcounts(self, joined, columns, output = None):
        if output is None: output = joined.replace('.tsv', '_normalized.tsv')
        working_dir = '/'.join(joined.split('/')[:-1])
        info = pd.read_csv(joined, sep = '\t')
        info[columns] = info[columns].fillna(value=0)
        info[columns].to_csv(working_dir + '/to_normalize.tsv', sep = '\t', index = False)
        print('Normalizing {} on columns {}'.format(joined, ','.join(columns)))
        self.run_command('Rscript MOSCA/scripts/normalization.R --readcounts {0}/to_normalize.tsv --output {0}/normalization_factors.txt'.format(
                working_dir))
        factors = open(working_dir + '/normalization_factors.txt').read().split('\n')[:-1]      # there is always the \n as last element
        
        for i in range(len(columns)):
            info[columns[i] + '_normalized'] = info[columns[i]] * float(factors[i])
        return info
    
    def check_bowtie2_index(self, index_prefix):
        files = glob.glob(index_prefix + '*.bt2')
        if len(files) < 6:
            return False
        return True
    
    '''
    Input:
        reference: name of contigs file from assembly
        reads: list, [forward reads, reverse reads]
        output: basename of outputs
        threads: number of threads to use
        blast: name of BLAST file
    Output:
        Will generate a bowtie2 index named contigs.replace(.fasta,_index), 
        GFF annotation file named blast.replace(.blast,.gff)
        SAM alignment and READCOUNTS files named output + .sam and .readcounts
    '''
    def perform_alignment(self, reference, reads, basename, threads = 1, blast = None,
                          blast_unique_ids = True, attribute = 'gene_id'):
        
        if not self.check_bowtie2_index(reference.replace('.fasta', '_index')):
            print('INDEX files not found. Generating new ones')
            self.generate_mg_index(reference, reference.replace('.fasta', '_index'))
        else:
            print('INDEX was located at ' + reference.replace('.fasta', '_index'))
        self.align_reads(reads, reference.replace('.fasta', '_index'), basename + '.sam',
                         basename + '_bowtie2_report.txt', log = basename + '.log', 
                         threads = threads)
        
        if blast is None:
            if not os.path.isfile(reference.replace('.fasta', '.gff')):
                print('GFF file not found at ' + reference.replace('.fasta','.gff') + 
                      '. Generating a new one.')
                self.build_gff_from_contigs(reference, 
                                            reference.replace('.fasta','.gff'))
            else:
                print('GFF file was located at ' + reference.replace('.fasta','.gff'))
        else:
            if not os.path.isfile(blast.replace('.blast', '.gff')):
                print('GFF file not found at ' + blast.replace('.blast', '.gff') + 
                      '. Generating a new one.')
                if blast_unique_ids:
                    self.build_gff(blast, blast.replace('.blast', '.gff'))
            else:
                print('GFF file was located at ' + blast.replace('.blast', '.gff'))
        
        
        self.run_command('htseq-count -i {} {}.sam {}{}'.format(attribute, basename, 
            (reference.replace('.fasta','.gff') if blast is None else blast.replace('.blast', '.gff')), 
            ('' if blast is not None else ' --stranded=no')), file = basename + '.readcounts')
    
    '''
    Input:
        df: pandas.DataFrame to manipulate
        column: column composed of lists from where to expand the dataframe
    Output:
        Returns the DataFrame expanded through one column by repeating all the
        values of the row where in 'column' there is a list with more than one 
        element
    '''    
    def expand_by_list_column(self, df, column = 'Pathway'):
        lens = [len(item) for item in df[column]]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary[column] = np.concatenate(df[column].values)
        return pd.DataFrame(dictionary) 
    
    def generate_mg_index(self, reference, index_prefix):
        self.run_command('bowtie2-build ' + reference + ' ' + index_prefix)
    
    def align_reads(self, reads, index_prefix, sam, report, log = None, threads = 6):  
        self.run_command('bowtie2 -x ' + index_prefix + ' -1 ' + reads[0] + 
                         ' -2 ' + reads[1] + ' -S ' + sam + ' -p ' + str(threads)
                         + ' 1> ' + report + ' 2> ' + log)
    
    def sort_alphanumeric(self, alphanumeric_list):
        return sorted(alphanumeric_list, key=lambda item: (int(item.partition(' ')[0])
                if item[0].isdigit() else float('inf'), item))
    
    def remove_files(self, files):
        for file in files:
            print('Deleting file', file)
            os.remove(file)
    
    def run_command(self, bashCommand, file = '', mode = 'w', sep = ' ', print_message = True, verbose = True):
        if print_message:
            print(bashCommand.replace(sep, ' '))
        if file == '':
                subprocess.run(bashCommand.split(sep), stdout=sys.stdout if verbose else None, 
                               check = True)
        else:
            with open(file, mode) as output_file:
                subprocess.run(bashCommand.split(sep), stdout=output_file)
    
    def run_pipe_command(self, bashCommand, output = '', mode = 'w', sep = ' ', print_message = True):
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
      
    def parse_blast(self, blast):
        result = pd.read_csv(blast, sep='\t', header = None)
        result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                          'evalue', 'bitscore']
        return result
    
    '''
    Input:
        file - str, name of blast file
        output - str, name of blast processed file to output
    Output:
        This function applies for blast files with multiple IDs for same protein.
        Only the first identification will be kept in a new file, to serve as 
        reference for GFF generation.
    '''
    def remove_multiple_ids_from_blast(self, file, output):
        print('Removing multiple identifications from ' + file + '. This may take a while...')
        file = self.parse_blast(file)
        output = pd.DataFrame(columns = file.columns)
        i = 0
        pbar = tqdm(total = len(file) + 1)
        while i < len(file):
            qseqid = file.iloc[i]['qseqid']
            output = output.append(file.iloc[i])
            i += 1
            while file.iloc[i]['qseqid'] == qseqid: 
                i += 1
            pbar.update(1)
        pbar.close()
            
        output.to_csv(output, sep = '\t', index = False)
    
    def build_gff(self, blast, output, assembler = 'metaspades'):
        gff = pd.DataFrame()
        diamond = self.parse_blast(blast)
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
        gff.to_csv(output, sep = '\t', index=False, header=False)
        
    def build_gff_from_contigs(self, contigs, output, assembler = None):
        gff = pd.DataFrame()
        contigs = self.parse_fasta(contigs)
        contigs = [[k,v] for k,v in contigs.items()]
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
        gff.to_csv(output, sep = '\t', index=False, header=False)
    
    def count_reads(self, file):
        lines = 0
        handler = gzip.open(file)
        for line in handler:
            lines += 1
        handler.close()
        return lines / 4
    	
    def parse_metaquast(self, file):
        result = dict()
        handler = open(file)
        handler.readline()
        flag = True
        while flag == True:
            line = handler.readline().rstrip('\n')
            print(line)
            if line == '':
                flag = False
            else:
                if line.startswith('#'):
                    result[line.split('\t')[0]] = line.split('\t')[1]
                else:
                    result[line.split('\t')[0]] = line.split('\t')[1]
        handler.close()
        return result
        
    def avaliate_annotation(self, ann, not_ann):
        annlines = len(open(ann).readlines())
        not_annlines = len(open(not_ann).readlines())
        return (annlines, not_annlines / 2, annlines / (not_annlines / 2))
    
    def merge_fq(self, file1, file2, output):
        self.run_command('bash MOSCA/scripts/merge-paired-reads.sh ' + file1 + 
                         ' ' + file2 + ' ' + output)
    
    def divide_fq(self, file, output1, output2):
        self.run_command('bash MOSCA/scripts/unmerge-paired-reads.sh ' + file + 
                         ' ' + output1 + ' ' + output2)
    def parse_fasta(self, file):
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
    
    def validate_arguments(self, parser):
        args = parser.parse_args()
        if args.files is None:
            if args.files is None:
                print('Must specify which files to use!')
            if args.output is None:
                print('Must specify which output directory to use!')
            parser.print_help()
            exit()
        args.output = args.output.rstrip('/')                                   # remove the automatic ending
        return args
    
    '''
    Input:
        experiment: list - filename(s) of files with reads
        type_of_data: 'mg' or 'mt'
        output: output directory
        sequencing_technology: 'paired' or 'single'
    Output:
        If data is interleaved (length of experiment is 1, and sequencing_technology
        is 'paired') it will divide reads to two files.
        Returns data, the files to perform the following steps in, and name, to
        be used as mg_name or mt_name
    '''  
    def process_argument_file(self, experiment, type_of_data, output,
                              sequencing_technology):
        data = experiment.split(',')
        
        if len(data) == 1 and sequencing_technology == 'paired':                  # if data is interleaved paired end, it will be split up to forward and reverse files
            name = data[0].split('/')[-1].split('.fastq')[0]
            (forward, reverse) = ('{}/Preprocess/{}{}'.format(output, name, fr) 
            for fr in ['_R1.fastq','_R2.fastq'])
            if os.path.isfile(forward):
                print('Split has already occurred for this file!')
            else:
                self.timed_message('Splitting reads from {} to {} and {}'.format(
                        data[0], forward, reverse))
                self.divide_fq(data[0], forward, reverse)
            data = [forward, reverse]
        elif len(data) == 1 and not sequencing_technology == 'paired':
            name = data[0].split('/')[-1].split('.fastq')[0]
        elif len(data) == 2 and sequencing_technology == 'paired':
            name = data[0].split('/')[-1].split('_R')[0]
        else:
            print("""BAD INPUT: {} reads must be given in one of three formats:
                    1. Paired-end reads - two files, R1 (forward) and R2(reverse)
                    2. Paired-end interleaved reads - one file
                    3. Single-end reads""".format(type_of_data.upper()))
            exit()
        return (data, name)

    def print_arguments(self, args):
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(vars(args))
        dict_args = vars(args)
        for key in dict_args.keys():
            key_name = key[0].upper() + key[1:].replace('_',' ')
            if key != 'files':
                print(key_name + ': ' + str(dict_args[key]))
            else:
                print(key_name)
                i = 0
                experiments = [experiment.split(':') for experiment in dict_args[key]]
                for experiment in experiments:
                    i += 1
                    print('Experiment' + str(i))
                    print('MG files: ' + experiment[0])
                    if len(experiment) > 1:
                        print('MT files: ' + experiment[1])
        
    '''
    input: a FASTA file of proteins identified by FGS, which might have some bases called as *
        a filename (temp) where to store the resulting file before moving it to the original file place
    output: a FASTA file where no protein sequence contains *
    '''
    def correct_fasta_file(self, fasta, temp = 'temp.fasta'):
        file = self.parse_fasta(fasta)
        handler = open(temp, 'w')
        for key,value in file.items():
            if value[0] == '*':             #FGS assigns * for the first position in a very small number of proteins
                value = value[1:]
            if '*' not in value:
                handler.write('>' + key + '\n' + value + '\n')
        os.rename(temp, fasta)
        
    '''
    Input:
        message: a message to be printed
    Output:
        will print the message with the time in human readable format
    '''
    def timed_message(self, message = None):
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ': ' + message)
            
    '''
    Input:
        task: Str, task that has finished ['preprocessing','assembly','annotation',
        'binning','expression']
        file: name of monitorization file
        task_output: name of folder/file where outputs are stored
    Output:
        will print the message with the time in human readable format
    '''
    def task_is_finished(self, task = None, file = None, task_output = None):
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + 
              ': {} has finished and results are available at {}'.format(
                task, task_output))
        with open(file, 'a') as f:
            if task in ['Metatranscriptomics analysis', 'Metaproteomics analysis']:
                task = 'Expression'
            f.write(task.lower() + '\n')
    
    '''
    Input:
        output_dir: str - the base project directory
        output_level: str - maximum, medium or minimum
    Output:
        All intermediate files of preprocessing will be removed if output_level < max,
        all files will be removed if output_level < med
    '''
    def remove_preprocessing_intermediates(self, output_dir, output_level):
        if output_level != 'maximum':
            for file in os.listdir(output_dir + '/FastQC'):
                if not any(s in file for s in ['quality_trimmed', '.html']):
                    file = '{}/FastQC/{}'.format(output_dir, file)
                    if os.path.isfile(file):
                        os.remove(file)
                    else:
                        shutil.rmtree(file)
            for file in (glob.glob(output_dir + '/SortMeRNA/*') + 
                         glob.glob(output_dir + '/Trimmomatic/*')):
                if not any(s in file for s in ['quality_trimmed', 'adapters.txt', 
                                           'quality_params.txt']):
                    os.remove(file)
            if output_level != 'medium':
                for file in (glob.glob(output_dir + '/FastQC/*.html') +
                             glob.glob(output_dir + '/Trimmomatic/*')):
                    os.remove(file)
                
    '''
    Input:
        output_dir: str - the base project directory
        output_level: str - maximum, medium or minimum
    Output:
        All intermediate files of assembly will be removed if output_level < max,
        all files will be removed if output_level < med
    '''
    def remove_assembly_intermediates(self, output_dir, output_level, samples):
        if output_level != 'maximum':
            for sample in samples:
                for file in os.listdir('{}/Assembly/{}'.format(output_dir, sample)):
                    if not any(s in file for s in ['quality_control', 'contigs.fasta']):
                        file = '{}/Assembly/{}/{}'.format(output_dir, sample, file)
                        if os.path.isfile(file):
                            os.remove(file)
                        else:
                            shutil.rmtree(file)
            if output_level != 'medium':
                for sample in samples:
                    os.remove('{}/Assembly/{}/contigs.fasta'.format(output_dir, sample))
                    shutil.rmtree('{}/Assembly/{}/quality_control'.format(output_dir, sample))

    '''
    Input:
        output_dir: str - the base project directory
        output_level: str - maximum, medium or minimum
    Output:
        All intermediate files of assembly will be removed if output_level < max,
        all files will be removed if output_level < med
    '''
    def remove_annotation_intermediates(self, output_dir, output_level, samples):
        if output_level != 'maximum':
            for sample in samples:
                for file in os.listdir('{}/Annotation/{}'.format(output_dir, sample)):
                    if not any(s in file for s in ['.faa', '.blast', 'cog', 'results']):
                        file = '{}/Annotation/{}/{}'.format(output_dir, sample, file)
                        os.remove(file)
            if output_level != 'medium':
                for sample in samples:
                    os.remove('{}/Annotation/{}/{}{}'.format(output_dir, sample, termination)
                    for termination in ['_fgs.faa', '_aligned.blast'])
                    
    '''
    Input:
        output_dir: str - the base project directory
        output_level: str - maximum, medium or minimum
    Output:
        All intermediate files of assembly will be removed if output_level < max,
        all files will be removed if output_level < med
    '''
    def remove_kegg_pathway_intermediates(self, output_dir, output_level, metabolic_map):
        if output_level != 'maximum':
            for file in glob.glob('{}/*{}*'.format(output_dir, metabolic_map)):
                os.remove(file)
    
    '''
    Input:
        filename: str - filename of FastQC report
    Output:
        returns dict{module:(value, pd.DataFrame)} with data from FastQC report
    '''      
    def parse_fastqc_report(self, filename):
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
                    partial_data = np.reshape(partial_data,(int(partial_data.size/len(labels)),len(labels)))
                    data[name] = (flag, pd.DataFrame(data = partial_data[1:,1:],
                                                    index = partial_data[1:,0],
                                                    columns = partial_data[0,1:]))
            i += 1
        return data
    
    '''
    Input:
        bashCommand: str - the command to retrieve the output from
        shell: bool - True if using some shell tool like grep or awk
    Output:
        Number of occurrences of character on file
    '''
    def count_on_file(self, expression, file, compressed = False):
        return int(subprocess.check_output("{} -c '{}' {}".format(
                'zgrep' if compressed else 'grep', expression, file), shell = True))
        
    '''
    Input:
        file: str - file to count lines of
    Output:
        Number of lines in file
    '''
    def count_lines(self, file):
        return int((subprocess.check_output("wc -l " + file, shell = True)).split()[0])
        
    '''
    Input:
        readcounts: str - file from htseq-count to normalize by contig size
        contigs: str - filename of contigs
    Output:
        The readcounts by contig in the file will be normalized by contig size
    '''
    def normalize_mg_readcounts_by_size(self, readcounts, contigs):
        self.run_pipe_command('head -n -5 {}'.format(readcounts), 
            file = readcounts.replace('.readcounts', '_no_tail.readcounts'))
        self.run_pipe_command("seqkit fx2tab {} | sort | awk '{{print $1\"\\t\"length($2)}}' | join - {} | awk '{{print $1\"\\t\"$3/$2}}'".format(
            contigs, readcounts.replace('.readcounts', '_no_tail.readcounts')), 
            file = readcounts.replace('.readcounts', '_normalized.readcounts'))
    
    '''
    Input:
        output: str - filename of output
        data: pd.DataFrame - data to write on several sheets
        lines: int - number of lines per sheet (Excel's maximum is 1048575)
        index: bool - write index or not
    Output:
        data will be outputed through several sheets
    '''
    def multi_sheet_excel(self, output, data, sheet_name = 'Sheet', 
                          lines = 1000000, index = False):
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
        i = 0; j = 1
        while i + lines < len(data):
            data.iloc[i:(i + lines)].to_excel(writer, sheet_name='{}({})'.format(sheet_name, str(j)), index = index)
            j += 1
        data.iloc[i:len(data)].to_excel(writer, sheet_name='{}({})'.format(sheet_name, str(j)), index = index)
        writer.save()