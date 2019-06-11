# -*- coding: utf-8 -*-
'''
General tools for support of MOSCA's functionalities

By João Sequeira

Jun 2017
'''

from io import StringIO
import pandas as pd
import subprocess, glob, re, os, gzip

class MoscaTools:
    
    '''
    Input: 
        relation: pandas.DataFrame from Annotater.join_reports
        origin_of_data: 'metagenomics', 'metatranscriptomics', 'compomics'
        name: name of column to add to relation
        readcounts: readcounts file from quantification with htseq-count
        blast: blast annotation file
        readcounts_has_tail: boolean, if file is htseq-count expression matrix with
        those last 5 lines of general information
        readcounts_has_last_ids: boolean, if file still uses the third portion of
        UniProt IDs (tr|2nd|3rd). I used to use those, I'm sorry
    Output: 
        'relation' df will receive additional column with abundance/expression 
        information.This column will be named 'name', if no 'name' input is given, 
        will be named as the folder containing the blast file if input, else as 
        the folder containing the readcounts file
    '''
    def define_abundance(self, relation, origin_of_data = 'metagenomics', name = None,
                         readcounts = None, blast = None, readcounts_has_tail = True, 
                         readcounts_has_last_ids = False):
        if name is None and not (readcounts is None and blast is None): 
            name = blast.split('/')[-2] if blast is not None else readcounts.split('/')[-2]    # name is assumed to be the same as the folder containing the file
        print('Adding info on sample ' + name)
        if readcounts is not None:                                              # necessary for 'metagenomics' and 'metatranscriptomics' (probably will be for all in the future)
            readcounts = pd.read_csv(readcounts, sep = '\t', header = None)  
            if readcounts_has_tail: readcounts = readcounts[:-5]                  # remove those last resume lines from htseq-count
            readcounts['Protein ID'] = ([ide.split('_')[0] for ide in readcounts[0]] 
                        if readcounts_has_last_ids else readcounts[0])
            readcounts.columns = ['geneid', name, 'Protein ID']
            del readcounts['geneid']
        if origin_of_data == 'metagenomics':                                    # quantifying 'metagenomics' makes use of quality control SAM file
            blast = self.parse_blast(blast)
            blast['contig'] = ['_'.join(ide.split('_')[:-3]) for ide in blast.qseqid]
            blast = pd.merge(blast, readcounts, left_on = 'contig', right_on = 'Protein ID')    # Protein ID are the contigs here. Don't judge me, my eyes are burning
            blast['Protein ID'] = [ide.split('|')[1] if ide != '*' else ide for ide in blast.sseqid]
            blast = blast.groupby('Protein ID')[name].sum().reset_index()[['Protein ID', name]]
            result = pd.merge(relation, blast, on = 'Protein ID', how = 'outer')
        elif origin_of_data == 'metatranscriptomics':
            result = pd.merge(relation, readcounts, on = 'Protein ID', how = 'outer')
        elif origin_of_data == 'compomics':
            pass
        result[name] = result[name].fillna(value = 0).astype(int)
        return result
    
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
        print('Normalizing ' + joined + ' on columns ' + ','.join(columns))
        self.run_command('Rscript MOSCA/tmm.R --readcounts ' + working_dir + '/to_normalize.tsv --output '
                         + working_dir + '/normalization_factors.txt')
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
    def perform_alignment(self, reference, reads, basename, threads = 1, blast = None):
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
                self.build_gff(blast, blast.replace('.blast', '.gff'))
            else:
                print('GFF file was located at ' + blast.replace('.blast', '.gff'))
        self.run_htseq_count(basename + '.sam', reference.replace('.fasta','.gff')
                            if blast is None else blast.replace('.blast', '.gff'),
                            basename + '.readcounts', 
                            stranded = False if blast is None else True)
    
    def generate_mg_index(self, reference, index_prefix):
        self.run_command('bowtie2-build ' + reference + ' ' + index_prefix)
    
    def align_reads(self, reads, index_prefix, sam, report, log = None, threads = 6):  
        self.run_command('bowtie2 -x ' + index_prefix + ' -1 ' + reads[0] + 
                         ' -2 ' + reads[1] + ' -S ' + sam + ' -p ' + str(threads)
                         + ' 1> ' + report + ' 2> ' + log)

    def run_htseq_count(self, sam, gff, output, attribute = 'gene_id', stranded = True):
        self.run_command('htseq-count -i ' + ' '.join([attribute, sam, gff]) + 
                         ('' if stranded else ' --stranded=no'), file = output)
    
    def sort_alphanumeric(self, alphanumeric_list):
        return sorted(alphanumeric_list, key=lambda item: (int(item.partition(' ')[0])
                if item[0].isdigit() else float('inf'), item))
    
    def remove_files(self, files):
        for file in files:
            print('Deleting file', file)
            os.remove(file)
    
    def run_command(self, bashCommand, file = '', mode = 'w', sep = ''):
        print(bashCommand)
        if file == '':
                subprocess.run(bashCommand.split(sep), stdout=subprocess.PIPE, check = True)       # was subprocess.Popen
        else:
            with open(file, mode) as output_file:
                subprocess.run(bashCommand.split(sep), stdout=output_file)           # was subprocess.Popen
        
    def parse_blast(self, blast):
        result = pd.read_csv(blast, sep='\t', header = None)
        result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                          'evalue', 'bitscore']
        return result
    
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
    	
    def meta_quast(self, fil,file,k):
        self.run_command('metaquast.py ' + file + ' -o Assembly/' + fil + '/' + k)
            
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
    
    def download_database(self, url):
        import os
        adress = url.replace('https','rsync')
        output = 'Databases/DIAMOND/' + url.split('genomes/')[1]
        if not os.path.isdir(output):
            os.makedirs(output)
        self.run_tool('rsync --copy-links --recursive --times --verbose ' + adress + ' ' + output)
        
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
        
    def uniprot_request(self, ids, output):
        import requests, time
         
        BASE_URL = 'http://www.uniprot.org/uniprot/'
         
        def http_get(url, params=None, stream=False):
            response = requests.get(url, params=params, stream=stream)
            return validate_response(response, stream)
         
        def validate_response(response, raw=False):
            if response.ok:
                if raw:
                    return response
                else:
                    return response.content
            else:
                response.raise_for_status()
         
        def http_post(url, params=None, json=None, headers=None, stream=False):
            response = requests.post(url, data=params, json=json,
                                     headers=headers, stream=stream)
            return validate_response(response, stream)
        
        print('tab')
        try:
            params = {
                'from':'ACC+ID',
                'format':'tab',
                'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                'columns':'id,ec,pathway,protein names'
            }
            response = http_post(BASE_URL, params=params)
            with open(output + '.tab', 'a') as f:
                f.write(response.decode('utf-8'))
            time.sleep(3)
        except:
            print(ids[0],'to',ids[-1],'tab failed')
            time.sleep(120)
            try:
                params = {
                    'from':'ACC+ID',
                    'format':'tab',
                    'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                    'columns':'id,ec,lineage(SUPERKINGDOM),lineage(PHYLUM),lineage(CLASS),lineage(ORDER),lineage(FAMILY),lineage(GENUS),lineage(SPECIES),pathway,protein names'
                }
                response = http_post(BASE_URL, params=params)
                with open(output + '.tab', 'a') as f:
                    f.write(response.decode('utf-8'))
                time.sleep(3)
            except:
                with open(output + '.log','a') as f: f.write(ids[0] + ' to ' + ids[-1] + ' tab failed again\n')
        
        for task in ['fasta','gff']:
            print(task)
            try:
                params = {
                    'from':'ACC+ID',
                    'format':task,
                    'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                }
                response = http_post(BASE_URL, params=params)
                with open(output + '.' + task, 'a') as f:
                    f.write(response.decode('utf-8'))
                time.sleep(3)
            except:
                print(ids[0],'to',ids[-1],task,' failed')
                time.sleep(120)
                try:
                    params = {
                        'from':'ACC+ID',
                        'format':task,
                        'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                    }
                    response = http_post(BASE_URL, params=params);print(response)
                    with open(output + '.' + task, 'a') as f:
                        f.write(response.decode('utf-8'))
                    time.sleep(3)
                except:
                    with open(output + '.log','a') as f: f.write(ids[0] + ' to ' + ids[-1] + ' ' + task + ' failed again\n')

    def chunky(self, ids, output, chunk):
        i = 0
        j = 0
        while j < len(ids):
            if i + chunk > len(ids):
                j = len(ids)
            else:
                j = i + chunk
            self.uniprot_request(ids[i:j], output)
            i = i + chunk
    
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
        if args.files is None or args.output is None:
            if args.files is None:
                print('Must specify which files to use!')
            if args.output is None:
                print('Must specify which output directory to use!')
            parser.print_help()
            exit()
        else:
            args.output = args.output.rstrip('/')                               # remove the automatic ending
            return args

    def print_arguments(self, args):
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(vars(args))
        print()
        dict_args = vars(args)
        print(dict_args)
        for key in dict_args.keys():
            key_name = key[0].upper() + key[1:].replace('_',' ')
            if key is not 'files':
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
        
def compare_result_with_reality(result, reality):       #result is uniprot tab, reality is tsv taxonomy1\ttaxonomy2\t...\tpercentage
    import pandas as pd    
    print('tis a mistery')
    resultdf = pd.read_csv(result, sep = '\t', header = None)
    resultdf.columns = ['id','ec','lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)','lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','pathway','protein names']
    resultdf = resultdf.groupby(resultdf.columns[2:9].tolist()).size().reset_index().rename(columns={0:'abundance'})
    realitydf = pd.read_csv(reality, sep = '\t', header = None)
    realitydf.columns = ['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                'lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','abundance']
    realitydf['lineage(SPECIES)'] = realitydf['lineage(GENUS)'] + ' ' + realitydf['lineage(SPECIES)']
    total = float(resultdf['abundance'].sum())
    print('Total is ' + str(total))
    resultdf['abundance'] /= (total/100)      #set as percentage
    totalreal = float(realitydf['abundance'].sum())
    realitydf['abundance'] /= (totalreal/100)
    #ammendments on the result because of some uniprot facts
    resultdf.loc[resultdf['Taxonomic lineage (FAMILY)'].str.contains('Chloroflexaceae'),'Taxonomic lineage (FAMILY)'] = 'Chloroflexaceae'
    resultdf.loc[resultdf['Taxonomic lineage (SPECIES)'].str.contains('Methanosarcina mazei'),'Taxonomic lineage (SPECIES)'] = 'Methanosarcina mazei'
    realitydf.loc[realitydf['lineage(GENUS)'].str.contains('Methanosaeta'),'lineage(GENUS)'] = 'Methanothrix'
    
    for value in realitydf['lineage(SPECIES)'].get_values():
        print(value)
        previous_columns = ['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                            'lineage(FAMILY)','lineage(GENUS)']
        print(resultdf[resultdf['lineage(SPECIES)'] == value and resultdf[previous_columns] != realitydf[previous_columns]])
    '''
    for level in ['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                'lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)']:
        print('Level: ' + level)
        factors = list(realitydf[level].unique())
        value = 0
        others = 100
        for factor in factors:
            print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['abundance'].sum())) + '%')
            others -= float(resultdf[resultdf[level] == factor]['abundance'].sum())
            value += (float(resultdf[resultdf[level] == factor]['abundance'].sum()) - abs(float(resultdf[resultdf[level] == factor]['abundance'].sum()) - realitydf.loc[realitydf[level] == factor]['abundance'].sum()))
        if level == 'lineage(SUPERKINGDOM)':
            for factor in ['Eukaryota','Viruses']:
                print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['abundance'].sum())) + '%')
                others -= float(resultdf[resultdf[level] == factor]['abundance'].sum())
        print('Others:' + str(others) + '%')
        print('Correctness: ' + str(value) + '%')
    '''

def compare_result_with_reality_coverage_tax(result, reality):       #result is uniprot tab, reality is tsv taxonomy1\ttaxonomy2\t...\tpercentage
    import pandas as pd    
    resultdf = pd.read_csv(result, sep = '\t')
    resultdf = resultdf[['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                         'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                         'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                         'Taxonomic lineage (SPECIES)','coverage']]
    resultdf = resultdf.groupby(resultdf.columns.tolist()[:-1])['coverage'].sum().reset_index()
    realitydf = pd.read_csv(reality, sep = '\t', header = None)
    realitydf.columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                         'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                         'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                         'Taxonomic lineage (SPECIES)','abundance']
    total = float(resultdf['coverage'].sum())
    print('Total is ' + str(total))
    resultdf['coverage'] /= (total/100)      #set as percentage
    totalreal = float(realitydf['abundance'].sum())
    realitydf['abundance'] /= (totalreal/100)
    #ammendments on the result because of some uniprot facts
    resultdf.loc[resultdf['Taxonomic lineage (FAMILY)'].str.contains('Chloroflexaceae'),'Taxonomic lineage (FAMILY)'] = 'Chloroflexaceae'
    resultdf.loc[resultdf['Taxonomic lineage (SPECIES)'].str.contains('Methanosarcina mazei'),'Taxonomic lineage (SPECIES)'] = 'Methanosarcina mazei'
    for level in ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)']:
        print('Level: ' + level)
        factors = list(realitydf[level].unique())
        value = 0
        others = 100
        for factor in factors:
            print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['coverage'].sum())) + '%')
            others -= float(resultdf[resultdf[level] == factor]['coverage'].sum())
            value += (float(resultdf[resultdf[level] == factor]['coverage'].sum()) - abs(float(resultdf[resultdf[level] == factor]['coverage'].sum()) - realitydf.loc[realitydf[level] == factor]['abundance'].sum()))
        if level == 'Taxonomic lineage (SUPERKINGDOM)':
            for factor in ['Eukaryota','Viruses']:
                print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['coverage'].sum())) + '%')
                others -= float(resultdf[resultdf[level] == factor]['coverage'].sum())
        print('Others:' + str(others) + '%')
        print('Correctness: ' + str(value) + '%')
        
def taxonomy_of_others(result, reality, output):
    import pandas as pd    
    resultdf = pd.read_csv(result, sep = '\t')
    taxcolumns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)','coverage']
    resultdf = resultdf[taxcolumns]
    resultdf = resultdf.groupby(resultdf.columns.tolist()[:-1])['coverage'].sum().reset_index()
    realitydf = pd.read_csv(reality, sep = '\t', header = None)
    realitydf.columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)',
                         'Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                         'Taxonomic lineage (SPECIES)','abundance']
    realitydf['Taxonomic lineage (SPECIES)'] = realitydf['Taxonomic lineage (GENUS)'] + ' ' + realitydf['Taxonomic lineage (SPECIES)']
    #ammendments on the result because of some uniprot facts
    resultdf.loc[resultdf['Taxonomic lineage (FAMILY)'].str.contains('Chloroflexaceae'),'Taxonomic lineage (FAMILY)'] = 'Chloroflexaceae'
    resultdf.loc[resultdf['Taxonomic lineage (SPECIES)'].str.contains('Methanosarcina mazei'),'Taxonomic lineage (SPECIES)'] = 'Methanosarcina mazei'
    resultdf.loc[resultdf['Taxonomic lineage (GENUS)'].str.contains('Methanothrix'),'Taxonomic lineage (GENUS)'] = 'Methanosaeta'
    simulatedspecies = realitydf['Taxonomic lineage (SPECIES)'].tolist()
    othersdf = resultdf[~resultdf['Taxonomic lineage (SPECIES)'].isin(simulatedspecies)]
    othersdf.to_excel(output)
        
def split(pathway):
    pathway = pathway.split('. ')
    return [path for path in pathway if path != '']

def using_repeat(df):
    import numpy as np
    import pandas as pd
    lens = [len(item) for item in df['Pathway']]
    dictionary = dict()
    for column in df.columns:
        dictionary[column] = np.repeat(df[column].values,lens)
    dictionary["Pathway"] = np.concatenate(df['Pathway'].values)
    return pd.DataFrame(dictionary)

def organize_pathway(result):
    import pandas as pd    
    import numpy as np
    #resultdf = pd.read_csv(result, sep = '\t')
    resultdf = result[result.Pathway.notnull()]
    resultdf.Pathway = resultdf.Pathway.apply(split)
    resultdf = using_repeat(resultdf)
    pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) for path in resultdf.Pathway], index = resultdf.index)
    resultdf = resultdf[['coverage']]
    pathways.columns = ['superpathway','pathway','subpathway']
    resultdf = pd.concat([pathways, resultdf], axis = 1)
    resultdf = resultdf.groupby(resultdf.columns.tolist()[:-1])['coverage'].sum().reset_index()
    return resultdf     #.set_index(resultdf.columns.tolist()[:-1])

def get_ids_from_ncbi_gff(gff):
    'RefSeq:'
    lines = [line.rstrip('\n') for line in open(gff).readlines() if 'CDS' in line]
    df = pd.read_csv(StringIO('\n'.join(lines)), sep = '\t', header=None)
    df.columns = ["seqid","source","type","start","end","score","strand","phase","attributes"]
    ids = list()
    for attribute in df.attributes:
        if 'RefSeq:' in attribute:
            ids.append(re.split(';|,', attribute.split('RefSeq:')[-1])[0])
        elif 'Genbank:' in attribute:
            ids.append(re.split(';|,', attribute.split('Genbank:')[-1])[0])
        else:
            ids.append(re.split(';|,', attribute.split('Dbxref=Genbank:')[-1])[0])
    return ids

if __name__ == '__main__':
    mtools = MoscaTools()
    
    for n in ['1','2','3','4']:
        reads = ['MOSCAfinal/Preprocess/Trimmomatic/quality_trimmed_4478-R' +
                 n + '-1-MiSeqKapa_' + fr + '_paired.fq' for fr in ['forward','reverse']]
        
        mtools.perform_alignment('MOSCAfinal/Assembly/joined/contigs.fasta', 
                                 reads, 'MOSCAfinal/Metatranscriptomics/mt' + n, 
                                 threads = 10)
        
'''
for name in ['OLDES6_S4_L001','PAL6_S2_L001']:
    pbar = ProgressBar()
    fasta = pd.DataFrame.from_dict(parse_fasta(name + '_database.fasta'), orient = 'index')
    fasta['qseqid'] = fasta.index
    fasta.reset_index(inplace = True)
    data = pd.merge(fasta, blast[['qseqid','id']], on = 'qseqid', how = 'inner')
    with open(name + '_uniprotids_database.fasta','w') as f:
        for i in pbar(range(len(data))):
            f.write('>' + data.iloc[i]['id'] + '\n' +  data.iloc[i][0] + '\n')
'''