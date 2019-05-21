# -*- coding: utf-8 -*-
"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

from diamond import DIAMOND
from mosca_tools import MoscaTools
from uniprot_mapping import UniprotMapping
from progressbar import ProgressBar
from io import StringIO
from tqdm import tqdm
import pandas as pd
import numpy as np
import time, os, shutil, glob

mtools = MoscaTools()
upmap = UniprotMapping()

class Annotater:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    '''
    Input:
        file: name of input file to perform gene calling on
        output: basename of output files
        assembled: True if input is contigs, False if it are reads
        error_model: quality model to consider when input are reads
    Output:
        FragGeneScan output files will be produced with basename 'output + /fgs'
        If input is FASTQ reads (if assembled == False) a FASTA version will
        be produced in the same folder with the same name as the FASTQ file
        (but with .fasta instead of .fastq)
    '''
    def gene_calling(self, file, output, assembled = True, error_model = 'illumina_10'):
        bashCommand = 'run_FragGeneScan.pl -genome='
        if assembled:
            bashCommand += file + ' -out=' + output + '/fgs -complete=1 -train=./complete'
        else:
            mtools.fastq2fasta(file, file.replace('fastq', 'fasta'))            # fraggenescan only accepts FASTA input
            bashCommand += (file.replace('fastq', 'fasta') + ' -out=' + output +
                            '/fgs -complete=0 -train=./' + error_model)
        bashCommand += ' -thread=' + self.threads
        mtools.run_command(bashCommand)
    
    def annotation(self):
        diamond = DIAMOND(threads = self.threads,
                          db = self.db,
                          out = self.out_dir + '/Annotation/' + self.name + '/aligned.blast',
                          query = self.out_dir + '/Annotation/' + self.name + '/fgs.faa',
                          un = self.out_dir + '/Annotation/' + self.name + '/unaligned.fasta',
                          unal = '1',
                          max_target_seqs = '1')
        
        if self.db[-6:] == '.fasta':
            print('FASTA database was inputed')
            if not os.path.isfile(self.db.replace('fasta','dmnd')):
                print('DMND database not found. Generating a new one')
                diamond.set_database(self.db, self.db.replace('fasta','dmnd'))
            else:
                print('DMND database was found. Using it')
        elif self.db[-5:] != '.dmnd':
            print('Database must either be a FASTA (.fasta) or a DMND (.dmnd) file')
            exit()
        diamond.db = self.db.split('.dmnd')[0] if '.dmnd' in self.db else self.db.split('.fasta')[0]
        
        diamond.run()
    
    def uniprot_request(self, ids, original_database = 'ACC+ID', database_destination = '',
                        output_format = 'tab', columns = None, databases = None):
        import requests
         
        BASE_URL = 'https://www.uniprot.org/uniprot/'
         
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
        
        params = {
            'from':original_database,
            'format':output_format,
            'query':'+OR+'.join(['accession:' + acc for acc in ids]),
            'columns':upmap.string4mapping(columns = columns, databases = databases)
        }
        if database_destination != '' or original_database != 'ACC+ID':
            params['to'] = 'ACC'
        try:
            return http_post(BASE_URL, params=params).decode('utf-8')
        except:
            try:
                return http_post(BASE_URL, params=params).decode('utf-8')
            except:
                return ''
            
    def get_uniprot_information(self, ids, original_database = 'ACC+ID', output_format = 'tab',
                                database_destination = '', chunk = 1000, sleep = 30,
                                columns = None, databases = None):
        pbar = ProgressBar()
        print('Retrieving UniProt information from ' + str(len(ids)) + ' IDs.')
        if output_format == 'tab':
            result = pd.DataFrame()
            for i in pbar(range(0, len(ids), chunk)):
                j = i + chunk if i + chunk < len(ids) else len(ids)
                data = self.uniprot_request(ids[i:j], original_database, database_destination)
                time.sleep(sleep)
                if len(data) > 0:
                    uniprot_info = pd.read_csv(StringIO(data), sep = '\t')
                    result = pd.concat([result, uniprot_info])
        elif output_format == 'fasta':
            result = str()
            for i in pbar(range(0, len(ids), chunk)):
                j = i + chunk if i + chunk < len(ids) else len(ids)
                data = self.uniprot_request(ids[i:j], original_database, 
                            database_destination, output_format = output_format, 
                            columns = columns, databases = databases)
                time.sleep(sleep)
                if len(data) > 0:
                    result += data
        return result
    
    def recursive_uniprot_fasta(self, output, fasta = None, blast = None, max_iter = 5):
        if fasta is not None:
            fasta = mtools.parse_fasta(fasta)
            all_ids = list(set(fasta.keys()))
        elif blast is not None:
            blast = mtools.parse_blast(blast)
            all_ids = list(set([ide.split('|')[1] if ide != '*' else ide 
                                for ide in blast.sseqid]))
        i = 0
        ids_done = ([ide.split('|')[1] for ide in mtools.parse_fasta(output).keys()]
                    if os.path.isfile(output) else list())
        while len(ids_done) < len(all_ids) and i < max_iter:
            print('Checking which IDs are missing information.')
            pbar = ProgressBar()
            ids_missing = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
            print('Information already gathered for ' + str(len(ids_done)) + 
                  ' ids. Still missing for ' + str(len(ids_missing)) + '.')
            uniprotinfo = self.get_uniprot_information(ids_missing, 
                                                       output_format = 'fasta')
            with open(output, 'a') as file:
                file.write(uniprotinfo)
            ids_done = [ide.split('|')[1] for ide in mtools.parse_fasta(output).keys()]
            i += 1
        if len(ids_done) == len(all_ids):
            print('Results for all IDs are available at ' + output)
        else:
            ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
            handler = open(ids_unmapped_output, 'w')
            handler.write('\n'.join(ids_missing))
            print(str(i) + ' iterations were made. Results related to ' + str(len(ids_missing)) + 
                  ' IDs were not obtained. IDs with missing information are available' +
                  ' at ' + ids_unmapped_output + ' and information obtained is available' +
                  ' at ' + output)
        
    def recursive_uniprot_information(self, blast, output, max_iter = 5):
        if os.path.isfile(output):
            result = pd.read_csv(output, sep = '\t').drop_duplicates()
            ids_done = set(list(result['Entry']))
        else:
            print(output + ' not found.')
            ids_done = list()
            result = pd.DataFrame()
        all_ids = set([ide.split('|')[1] for ide in DIAMOND(out = blast).parse_result()['sseqid'] if ide != '*'])
        i = 0
        ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
        print('Checking which IDs are missing information.')
        pbar = ProgressBar()
        ids_missing = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
        print('ids missing:'+str(len(ids_missing)))
        print('ids done:'+str(len(ids_done)))
        print('ids all:'+str(len(all_ids)))
        
        while len(ids_missing) > 0 and i < max_iter:
            try:
                print('Information already gathered for ' + str(len(ids_done)) + 
                      ' ids. Still missing for ' + str(len(ids_missing)) + '.')
                uniprotinfo = self.get_uniprot_information(ids_missing)
                result = pd.concat([result, uniprotinfo])
                ids_done = set(list(result['Entry']))
                print('Checking which IDs are missing information.')
                pbar = ProgressBar()
                ids_missing = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
                i = 0
            except:
                print('Failed to retrieve information for some IDs. Retrying request.')
                i += 1
            
        result.to_csv(output, sep = '\t', index = False)
        
        if len(ids_missing) == 0:
            print('Results for all IDs are available at ' + output)
        else:
            handler = open(ids_unmapped_output, 'w')
            handler.write('\n'.join(ids_missing))
            print('Maximum iterations were made. Results related to ' + str(len(ids_missing)) + 
                  ' IDs were not obtained. IDs with missing information are available' +
                  ' at ' + ids_unmapped_output + ' and information obtained is available' +
                  ' at ' + output)
        
    '''
    Input:
        pathway: a String row from UniProt's 'Pathway' column
    Output:
        returns List of pathways included in that row
    '''
    def split(self, pathway):
        pathway = pathway.split('. ')
        return [path for path in pathway if path != '']
    
    '''
    Input:
        ec: a String row from UniProt's 'EC number' column
    Output:
        returns List of EC numbers included in that row
    '''
    def split_ec(self, ec):
        ec = ec.split('; ')
        return [ec_number for ec_number in ec if ec_number != '']
    
    '''
    Input:
        pathway: a String row from UniProt's 'Pathway' column
    Output:
        Reeives a pd.DataFrame object and breaks the list elements in the 'Pathway'
        column, multiplicating the rows with several elements of 'Pathway' into
        individual 'Pathway' elements, each with its row
    '''    
    def using_repeat(self, df):
        import numpy as np
        import pandas as pd
        lens = [len(item) for item in df['Pathway']]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary["Pathway"] = np.concatenate(df['Pathway'].values)
        return pd.DataFrame(dictionary) 
    
    '''
    Input:
        uniprotinfo: information from UniProt ID mapping
        blast: blast file from DIAMOND annotation
        output: basename for EXCEL files to output
    Output:
        Two EXCEL files formated for Krona plotting with taxonomic and functional 
        information.
        This function is very useful if wanting to use UniProt 'Pathway' information,
        as it uses the three previous functions to organize the information from
        that column into the three functional levels of UniProt Pathways.
        This function was very cool
    '''
    def uniprotinfo_to_excel(self, uniprotinfo, blast, output):
        blast = DIAMOND(out = blast).parse_result()
        uniprotdf = pd.read_csv(uniprotinfo, sep = '\t', index_col = 0).drop_duplicates()
        pbar = ProgressBar()
        blast['Coverage'] = [float(ide.split('_')[5]) for ide in pbar(blast.qseqid)]
        pbar = ProgressBar()
        blast['ID'] = [ide.split('|')[-1] for ide in pbar(blast.sseqid)]
        uniprotdf = pd.merge(uniprotdf, blast, left_index = True, right_on = 'ID')
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)','Coverage']
        taxdf = uniprotdf[tax_columns].groupby(tax_columns[:-1])['Coverage'].sum().reset_index()
        taxdf.to_excel(output + '_taxonomic_krona.xlsx', index = False)
        print('Saved taxonomy')
        func_columns = ['Pathway','Protein names','EC number']
        funcdf = uniprotdf[uniprotdf.Pathway.notnull()][func_columns + ['Coverage']]
        funcdf.Pathway = funcdf.Pathway.apply(self.split)
        funcdf = self.using_repeat(funcdf)
        pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) 
                                    for path in funcdf.Pathway], index = funcdf.index)
        pathways.columns = ['Superpathway','Pathway','Subpathway']
        del funcdf['Pathway']
        funcdf = pd.concat([pathways, funcdf], axis = 1)
        funcdf = funcdf[['Superpathway','Pathway','Subpathway','EC number',
                         'Protein names','Coverage']]
        funcdf = funcdf.groupby(funcdf.columns.tolist()[:-1])['Coverage'].sum().reset_index()
        funcdf.to_excel(output + '_functional_krona.xlsx', index = False)
        print('Saved pathways')
           
    def info_with_coverage_metaspades(self, blast, output):
        blast = DIAMOND(out = blast).parse_result()
        coverage = pd.Series([float(ide.split('_')[5]) for ide in blast.qseqid])
        blast['coverage'] = coverage.values
        coverage_values = list(set(blast.coverage))
        result = pd.DataFrame()
        for value in coverage_values:
            print(value)
            ids = [ide.split('|')[-1] for ide in blast[blast['coverage'] == value]['sseqid']]
            print(len(ids))
            try:
                part = self.get_uniprot_information(ids)
                part['coverage'] = pd.Series([value for i in range(len(part))]).values
                result = pd.concat([result, part])
            except:
                try:
                    part = self.get_uniprot_information(ids)
                    part['coverage'] = pd.Series([value for i in range(len(part))]).values
                    result = pd.concat([result, part])
                except:
                    result.to_csv(output, sep = '\t', index = False)
                    with open('errors.log', 'a') as f:
                        f.write(str(value) + ' failed')
                    continue
        result.to_csv(output, sep = '\t', index = False)
        
    def contigs_readcounts(self, contigs, output, read1, read2):
        self.run_alignment([[read1, read2]], contigs, output)
        
    def info_with_coverage_megahit(self, blast, readcounts, output, output_others = None):
        blast = DIAMOND(out = blast).parse_result()
        blast.qseqid = ['_'.join(qseqide.split('_')[:2]) for qseqide in blast.qseqid]
        readcountsdf = pd.read_csv(readcounts, sep = '\t', header = None)
        readcountsdf.columns = ['contig', 'coverage']
        pbar = ProgressBar()
        for value in pbar(readcountsdf.contig.get_values().tolist()):
            blast.loc[blast.qseqid == value,'coverage'] = float(readcountsdf[readcountsdf.contig==value]['coverage'])
        coverage_values = list(set(blast.coverage))
        result = pd.DataFrame()
        blast.to_csv('blast.blast', sep = '\t', index = False)
        for value in coverage_values:
            print(value)
            ids = [ide.split('|')[-1] for ide in blast[blast['coverage'] == value]['sseqid']]
            print(len(ids))
            try:
                part = self.get_uniprot_information(ids)
                part['coverage'] = pd.Series([value for i in range(len(part))]).values
                result = pd.concat([result, part])
            except:
                try:
                    part = self.get_uniprot_information(ids)
                    part['coverage'] = pd.Series([value for i in range(len(part))]).values
                    result = pd.concat([result, part])
                except:
                    result.to_csv(output, sep = '\t', index = False)
                    with open('errors.log', 'a') as f:
                        f.write(str(value) + ' failed')
                    continue
        result.to_csv(output, sep = '\t', index = False)

    def ni_proteins(self, fasta, blast, output, ni_proteins = True):
        proteins = mtools.parse_fasta(fasta)
        blast = DIAMOND(out = blast).parse_result()
        blast.index = blast.qseqid
        ids = list(blast[blast.sseqid == '*'].index) if ni_proteins == True else list(blast[blast.sseqid != '*'].index)
        handler = open(output, 'w')
        for ide in ids:
            if '*' not in proteins[ide]:
                handler.write('>' + ide + '\n' + proteins[ide] + '\n')
        handler.close()
        
    def parse_interproscan_output(self, file):
        df = pd.read_csv(file, sep = '\t', header = None)
        df.columns = ['Protein Accession', 'Sequence MD5 digest', 'Sequence Length', 
                      'Analysis', 'Signature Accession', 'Signature Description', 
                      'Start location', 'Stop location', 'Score', 'Status', 'Date', 
                      'InterPro annotations - accession', 'InterPro annotations - description', 
                      'GO annotations', 'Pathways annotations']
        return df
        
    def correct_interproscan_file(self, file):
        lines = open(file).readlines()
        lines = [line.rstrip('\n') + '\t' * (14-line.count('\t')) + '\n' for line in lines]
        handler = open(file, 'w')
        handler.write(''.join(lines))
        handler.close()
        
    '''
    Input: 
        blast: name of an annotated blast file, probably with a lot of not identified (*) proteins
        interproscan: name of an InterProScan result file with hopefully CDD annotations
        output: name of file to output
    Output: 
        an annotated blast file where, if the respectively protein had identifications in
        CDD database, now in the place of proteins names, it has the CDD IDs (may result in additional
        lines for such protein, since several domains might be identified) (output)
    '''        
    def blast2cddblast(self, blast, interproscan, output, correct_interproscan = True):
        blastdf = DIAMOND(out = blast).parse_result()
        if correct_interproscan: self.correct_interproscan_file(interproscan)
        interproscandf = self.parse_interproscan_output(interproscan)
        interproscandf = interproscandf[interproscandf['Analysis'] == 'CDD']
        cddblastdf = pd.DataFrame(columns = blastdf.columns)
        print('Building BLAST file with CDD annotations...')
        pbar = ProgressBar()
        for i in pbar(range(len(interproscandf))):
            line = blastdf.iloc[i]
            line['sseqid'] = interproscandf.iloc[i]['Signature Accession']
            cddblastdf = cddblastdf.append(line)
        cddblastdf.to_csv(output, index = False, header = False, sep = '\t')
        
    '''
    Input: name of a fasta file of proteins to be annotated
        COG blast DB namebase from ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz
        name of output file
    Output: annotated file with CDD IDs
    '''
    def run_rpsblast(self, fasta, output, cog, threads = '0'):
        bashCommand = ('rpsblast -query ' + fasta + ' -db "' + cog + '" -out ' +
                       output + ' -outfmt 6 -num_threads ' + threads + 
                       ' -max_target_seqs 1')
        open('MOSCA/Databases/COG/command.bash', 'w').write(bashCommand + '\n') # subprocess was not handling well running this command, so an intermediate file is written with the command to run
        print(' '.join(bashCommand.split(',')))
        mtools.run_command('bash MOSCA/Databases/COG/command.bash')
        os.remove('MOSCA/Databases/COG/command.bash')
        
    '''
    Input: 
        name of blast output with CDD IDs
        name of cddid summary file from ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
        name of fun file available at ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/fun.txt
        name of whog file available at ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG/whog
        name of cdd2cog script
        output folder where to store the resuls folder
    Output: results folder 
    '''
    def annotate_cogs(self, blast, output, cddid, fun, whog):
        mtools.run_command('perl MOSCA/scripts/cdd2cog.pl -r ' + blast + ' -c ' + 
                           cddid + ' -f ' + fun + ' -w ' + whog)
        if os.path.isdir(output + '/results'):
            shutil.rmtree(output + '/results')
        os.rename('results', output + '/results')
        
    '''
    Input: 
        cogblast: the output from cdd2go, a blast file with CDD and COG annotations
        fun: the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
    Output: 
        returns pandas.DataFrame with the functional categories intrisic levels 
        reorganized into corresponding columns
    '''      
    def organize_cdd_blast(self, cogblast, fun = 'MOSCA/Databases/COG/fun.txt'):
        cogblast = self.parse_cogblast(cogblast)
        cogblast = cogblast[cogblast['functional categories'].notnull()]
        cog_relation = self.parse_fun(fun)
        data = [cog_relation[functional_category] for functional_category in cogblast['functional categories']]
        result = pd.DataFrame(data)
        result.columns = ['COG general functional category','COG functional category']
        result = pd.concat([result[['COG general functional category','COG functional category']], 
                            cogblast[['COG protein description','cog','qseqid']]], axis = 1)
        return result
    
    '''
    Input: 
        cogblast: the output from cdd2go, a blast file with CDD and COG annotations
        output: filename of EXCEL file to write
    Output: 
        an EXCEL file with the COG identifications counted for krona plotting will
        be written
    '''      
    def write_cogblast(self, cogblast, output):
        cogblast = self.organize_cdd_blast(cogblast)
        del cogblast['qseqid']                                                  # TODO - this is something that should be reworked in the future. self.organize_cdd_blast is called twice, and while here the qseqid is not needed, it is needed in the self.join_reports call of self.global_information
        cogblast = cogblast.groupby(cogblast.columns.tolist()).size().reset_index().rename(columns={0:'count'})
        cogblast.to_excel(output, index = False)
        
    '''
    Input: name of cddblast to parse
    Output: pandas.DataFrame object
    '''      
    def parse_cogblast(self, cogblast):
        cogblast = pd.read_csv(cogblast, header=None, skiprows = 1, sep = '\t', low_memory=False)
        cogblast = cogblast[list(range(0,14))+[18]]                                             #several extra columns are produced because of bad formatting
        cogblast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                   'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'cog',
                   'functional categories', 'COG protein description']
        return cogblast
        
    '''
    Input: the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
    Output: a dictionary in the form {COG category (letter): (COG category (name), COG supercategory)}
    '''      
    def parse_fun(self, fun):
        lines = open(fun).readlines()
        result = dict()
        supercategory = lines[0]
        i = 1
        while i < len(lines):
            line = lines[i].rstrip('\n')
            if '[' in line:
                letter = line.split('[')[-1].split(']')[0]
                name = line.split('] ')[-1]
                result[letter] = [supercategory.rstrip('\n'), name]
            else:
                supercategory = line
            i += 1
        return result
    
    '''
    Input: 
        blast: name of the file from the DIAMOND annotation
        uniprotinfo: name of the uniprot information file
        cog_blast: name of the COG annotation file from self.annotate_cogs
        fun: name of the fun.txt file
        split_pathways: boolean, if MOSCA should split the information from the
        Pathway column provided by UniProt mapping
        readcounts_matrix: name of the file containing protein expression
    Output: 
        pandas.DataFrame with integration of information from UniProt and COGs
    '''
    def join_reports(self, blast, uniprotinfo, cog_blast, fun, split_pathways = False):
        result = mtools.parse_blast(blast)
        result.index = result.qseqid
        result = result[result.sseqid != '*']
        result['Entry'] = [ide.split('|')[1] for ide in result.sseqid]
        uniprotinfo = pd.read_csv(uniprotinfo, sep= '\t', low_memory=False).drop_duplicates()
        if split_pathways == True:                                              # TODO - the reorganization of pathways incurs in extra lines for same IDs. Find workaround
            print('Reorganizing pathways information.')
            funcdf = uniprotinfo[uniprotinfo.Pathway.notnull()][['Entry','Pathway']]
            funcdf.Pathway = funcdf.Pathway.apply(self.split)
            funcdf = self.using_repeat(funcdf)
            pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) 
                                        for path in funcdf.Pathway], index = funcdf.index)
            pathways.columns = ['Superpathway','Pathway','Subpathway']
            del funcdf['Pathway']; del uniprotinfo['Pathway']
            funcdf = pd.concat([pathways, funcdf], axis = 1)
            uniprotinfo = pd.merge(uniprotinfo, funcdf, on = ['Entry'], how = 'outer')
        result = pd.merge(result, uniprotinfo, on = ['Entry'], how = 'outer')
        cog_blast = self.organize_cdd_blast(cog_blast, fun)
        result = pd.merge(result, cog_blast, on = ['qseqid'], how = 'outer')
        print('Defining best consensus COG for each UniProt ID.')
        cogs_df = pd.DataFrame()
        tqdm.pandas()
        cogs_df['cog'] = result.groupby('Entry')['cog'].progress_apply(lambda x:
            x.value_counts().index[0] if len(x.value_counts().index) > 0 else np.nan)
        cogs_relation = result[['COG general functional category','COG functional category',
                            'COG protein description','cog']]
        cogs_df['Entry'] = cogs_df.index
        cogs_relation = cogs_relation.drop_duplicates()
        cogs_df = cogs_df[cogs_df['cog'].notnull()]
        cogs_df = pd.merge(cogs_df, cogs_relation, on = 'cog', how = 'inner')
        result.drop(['COG general functional category', 'COG functional category',
              'COG protein description', 'cog'], axis=1, inplace=True)
        result = pd.merge(result, cogs_df, on = 'Entry', how = 'outer')
        result = result[['Entry','Taxonomic lineage (SUPERKINGDOM)',
              'Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)',
              'Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)',
              'Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)',
              'EC number', 'Ensembl transcript','Function [CC]', 'Gene names',
              'Gene ontology (GO)', 'Keywords','Pathway', 'Protein existence', 
              'Protein families', 'Protein names','Cross-reference (BRENDA)', 
              'Cross-reference (BioCyc)','Cross-reference (CDD)', 
              'Cross-reference (InterPro)','Cross-reference (KEGG)', 
              'Cross-reference (KO)','Cross-reference (Pfam)', 
              'Cross-reference (Reactome)','Cross-reference (RefSeq)', 
              'Cross-reference (UniPathway)','Cross-reference (eggNOG)',
              'COG general functional category', 'COG functional category',
              'COG protein description','cog']]
        result.columns = ['Protein ID'] + result.columns.tolist()[1:]
        result = result.drop_duplicates()
        return result
    
    def run(self):
        self.gene_calling(self.file, self.out_dir + '/Annotation/' + self.name, self.assembled)
        self.annotation()
    
    '''
    Input: 
        faa: FASTA file with protein sequences to be annotated
        output: name of folder where to output results
        cog: name of COG
        cddid: name of cddid.tbl file for COG analysis
        whog: name of whog file for COG analysis
        fun: name of fun.txt file for COG analysis
        databases: LIST, name of databases in format for COG analysis
        threads: STR, number of threads to use
    Output: 
        pandas.DataFrame with abundance and expression information
    '''
    def cog_annotation(self, faa, output, cddid = 'MOSCA/Databases/COG/cddid.tbl',
                       whog = 'MOSCA/Databases/COG/whog', fun = 'MOSCA/Databases/COG/fun.txt', 
                       databases = None, threads = '8'):
        self.create_split_cog_db('MOSCA/Databases/COG', 'MOSCA/Databases/COG/Cog', threads = threads)
        pns = glob.glob('./MOSCA/Databases/COG/Cog_' + threads + '_*.pn')
        databases = [pn.split('.pn')[0] for pn in pns]
        self.run_rpsblast(faa, output + '/cdd_aligned.blast', ' '.join(databases))
        if os.path.isdir(os.getcwd() + '/results'):                              # the cdd2cog tool does not overwrite, and fails if results directory already exists
            print('Eliminating ' + os.getcwd() + '/results')
            shutil.rmtree(os.getcwd() + '/results', ignore_errors=True)          # is not necessary when running the tool once, but better safe then sorry!
        
        self.annotate_cogs(output + '/cdd_aligned.blast', output, cddid, fun, whog)
        self.write_cogblast(output + '/results/rps-blast_cog.txt', output + '/cogs.xlsx')

    def set_to_uniprotID(self, fasta, aligned, output):
        pbar = ProgressBar()
        result = DIAMOND(out = aligned).parse_result()
        sequences = mtools.parse_fasta(fasta)
        print('Changing names of ' + fasta + '\nto identifications in ' + aligned + '\nand outputing to ' + output)
        with open(output,'w') as f:
            for key, value in pbar(sequences.items()):
                try:
                    f.write('>' + str(result[result.qseqid == key]['sseqid'].item()) + '\n' + value + '\n')
                except:
                    print(result[result.qseqid == key]['sseqid'])
                    
    def global_information(self):
        # Join reports
        if not os.path.isfile(self.out_dir + '/Annotation/fgs.faa'):
            mtools.run_command('cat ' + ' '.join(glob.glob(self.out_dir + '/Annotation/*/fgs.faa')),
                                             self.out_dir + '/Annotation/fgs.faa')
        if not os.path.isfile(self.out_dir + '/Annotation/aligned.blast'):
            mtools.run_command('cat ' + ' '.join(glob.glob(self.out_dir + '/Annotation/*/aligned.blast')), 
                           file = self.out_dir + '/Annotation/aligned.blast')
        
        # Retrieval of information from UniProt IDs
        self.recursive_uniprot_information(self.out_dir + '/Annotation/aligned.blast', 
                                           self.out_dir + '/Annotation/uniprot.info')
        
        # Functional annotation with COG database
        self.cog_annotation(self.out_dir + '/Annotation/fgs.faa', 
                            self.out_dir + '/Annotation', self.cddid, self.whog, 
                            self.fun, threads = self.threads)

        # Integration of all reports - BLAST, UNIPROTINFO, COG
        joined = self.join_reports(self.out_dir + '/Annotation/aligned.blast', 
                               self.out_dir + '/Annotation/uniprot.info', 
                               self.out_dir + '/Annotation/results/rps-blast_cog.txt', 
                               self.fun)

        blast_files = glob.glob(self.out_dir + '/Annotation/*/aligned.blast')
        for file in blast_files:
            mg_name = file.split('/')[-2]
            mtools.build_gff_from_contigs(self.out_dir + '/Assembly/' + mg_name + '/contigs.fasta', 
                    self.out_dir + '/Assembly/' + mg_name + '/quality_control/alignment.gff')
            mtools.run_htseq_count(self.out_dir + '/Assembly/' + mg_name + '/quality_control/alignment.sam', 
                                   self.out_dir + '/Assembly/' + mg_name + '/quality_control/alignment.gff',
                                   self.out_dir + '/Assembly/' + mg_name + '/quality_control/alignment.readcounts',
                                   stranded = False)
            joined = mtools.define_abundance(joined, 
                            readcounts = self.out_dir + '/Assembly/' + mg_name + 
                            '/quality_control/alignment.readcounts', blast = file)
        if os.path.isfile(self.out_dir + '/joined_information.xlsx'):
            os.remove(self.out_dir + '/joined_information.xlsx')
        joined.to_csv(self.out_dir + '/joined_information.tsv', index=False, sep='\t')
        print('joined was written to ' + self.out_dir + '/joined_information.tsv')
        
        writer = pd.ExcelWriter(self.out_dir + '/joined_information.xlsx', 
                                engine='xlsxwriter')
        i = 0
        j = 1
        while i + 1048575 < len(joined):
            joined.iloc[i:(i + 1048575)].to_excel(writer, sheet_name='Sheet ' + str(j), index = False)
            j += 1
        joined.iloc[i:len(joined)].to_excel(writer, sheet_name='Sheet ' + str(j), index = False)
        print('joined was written to ' + self.out_dir + '/joined_information.xlsx')
        return joined

    '''
    UniProt regularly updates its databases. As we are working with MG here, many
    identifications will sometimes pertain to ORFs or pseudogenes that have been 
    wrongly predicted to code for proteins. It may also happen that the original
    authors decided to withdraw their published sequences.
    Input:
        uniprotinfo: name of UniProt info file
    Output:
        Rows in UniProt info file that lack a species identification will be 
        updated to include that new information
    '''
    def info_from_no_species(self, uniprotinfo):
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        missing_uniprotinfo = uniprotinfo[uniprotinfo['Taxonomic lineage (SPECIES)'].isnull()]
        ids = list(set([ide for ide in missing_uniprotinfo['Entry']]))
        new_uniprotinfo = self.get_uniprot_information(ids)
        print(missing_uniprotinfo)
        for entry in new_uniprotinfo['Entry']:
            missing_uniprotinfo[missing_uniprotinfo['Entry'] == entry] = new_uniprotinfo[new_uniprotinfo['Entry'] == entry]
        print(missing_uniprotinfo)
        print(new_uniprotinfo)
        print(new_uniprotinfo[new_uniprotinfo['Taxonomic lineage (SPECIES)'].notnull()]['Taxonomic lineage (SPECIES)'])
        new_uniprotinfo.to_csv('test.tsv', sep = '\t', index = False)
    
    '''
    Input:
        smp_directory: foldername where the SMP files are. These files are
        obtained from ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
        output: basename for PN and databases
        threads: STR, number of threads that the workflow will use
        step: number of SMP files per database
    Output:
        threads - 1 databases will be outputed, each with a consecutive part of
        the list of SMP files available. These databases are formated for RPS-BLAST
        search
    '''
    def create_split_cog_db(self, smp_directory, output, threads = '6', step = None):
        dbs = (open('MOSCA/Databases/COG/databases.txt').read().split('\n') if
        os.path.isfile('MOSCA/Databases/COG/databases.txt') else list())
        if threads in dbs:
            print('Already built COG database for this number of threads')
        else:
            print('Generating COG databases for [' + threads + '] threads.')
            smp_list = glob.glob(smp_directory + '/COG*.smp')
            if step is None:
                step = round(len(smp_list) / float(threads))
            i = 0
            pn_files = list()
            output += '_' + threads + '_'
            while i + step < len(smp_list):
                pn_files.append(output + str(int(i / step)) + '.pn')
                open(output + str(int(i / step)) + '.pn', 'w').write('\n'.join(smp_list[i:i + step]))
                i += step
            open(output + str(int(i / step)) + '.pn', 'w').write('\n'.join(smp_list[i:len(smp_list)]))
            pn_files.append(output + str(int(i / step)) + '.pn')
            for file in pn_files:
                mtools.run_command('makeprofiledb -in ' + file + ' -title ' + file.split('.pn')[0] + 
                                   ' -out ' + file.split('.pn')[0])             # -title and -out options are defaulted as input file name to -in argument; -dbtype default is 'rps'
            open('MOSCA/Databases/COG/databases.txt','w').write('\n'.join(dbs + [threads]))
    
    def kegg_mapper(self, ids, output, step = 50):
        result = pd.DataFrame()
        ids_failed = list()
        for i in range(0, len(ids) - step, step):
            try:
                url = 'http://rest.kegg.jp/conv/genes/uniprot:' + '+uniprot:'.join(ids[i:i+step])
                bashCommand = 'wget ' + url + ' -O ' + output
                mtools.run_command(bashCommand)
                result = pd.concat([result, pd.read_csv(output, sep='\t', header=None)])
                print(str(round(i/len(ids) * 100),2) + '% done')
                print('Already gathered ' + str(len(result)) + ' ids.')
            except:
                ids_failed += ids[i:i+step]
                print('Mapping failed for some IDs.')
        return result, ids_failed

    
if __name__ == '__main__':
    '''
    ids = DIAMOND(out = 'MOSCAfinal/Annotation/joined/aligned.blast').parse_result()['sseqid']
    
    ids = [ide.split('|')[1] for ide in ids if ide != '*']
    '''
    '''
    print(os.getcwd())
    ui = pd.read_csv('uniprot.info',sep='\t')
    found_ids = ui['Entry'].tolist()
    
    all_ids = open('ids_missing.txt').readlines()[0].rstrip('\n').split(',')
    missing_ids = [ide for ide in all_ids if ide not in found_ids]
    print('Found IDs: ' + str(len(found_ids)))
    print('Missing IDs: ' + str(len(missing_ids)))
    '''
    '''
    annotater = Annotater()
    result = annotater.recursive_uniprot_information('MOSCAfinal/Annotation/aligned.blast',
                                                     'MOSCAfinal/Annotation/uniprot.info1')
    '''
    #result = pd.concat([result, ui])
    #result.to_csv('uniprot.info',sep='\t',index=False)
    
    '''
    mosca_dir = os.path.dirname(os.path.realpath(__file__))
    
    annotater = Annotater(out_dir = 'MGMP',
                      fun = mosca_dir + '/Databases/COG/fun.txt',
                      cog = mosca_dir + '/Databases/COG/Cog',
                      cddid = mosca_dir + '/Databases/COG/cddid.tbl',
                      whog = mosca_dir + '/Databases/COG/whog',
                      cdd2cog_executable = mosca_dir + '/cdd2cog.pl')
    
    annotater.global_information()
    '''
    
    annotater = Annotater()
    
    annotater.create_split_cog_db('/home/jsequeira/COGsmp', 'MOSCA/Databases/COG/split')

