# -*- coding: utf-8 -*-
"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

from diamond import DIAMOND
from mosca_tools import MoscaTools
from progressbar import ProgressBar
from io import StringIO
import pandas as pd
import numpy as np
import time, os

mtools = MoscaTools()

class Annotater:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    def gene_calling(self):
        bashCommand = 'perl ' + os.path.expanduser('~/FGS/run_FragGeneScan.pl') + ' -genome='
        if self.assembled == True:
            fgs_dict = {'megahit':'final.contigs.fa','metaspades':'contigs.fasta'}
            bashCommand += self.out_dir + '/Assembly/' + self.name + '/' + fgs_dict[self.assembler]
            bashCommand += ' -out=' + self.out_dir + '/Annotation/' + self.name + '/fgs -complete=1 -train=./complete'
        else:
            fastq2fasta_command = ("paste - - - - < " + self.out_dir + "/Preprocess/SortMeRNA/rejected.fastq.fastq"
                                   + "| cut -f 1,2 | sed 's/^@/>/' | tr \"\t" "\n\" > " + self.out_dir + 
                                   "/Preprocess/SortMeRNA/rejected.fasta")
            mtools.run_command(fastq2fasta_command)
            bashCommand += (self.out_dir + '/Preprocess/SortMeRNA/rejected.fasta -out=' + self.out_dir + 
            '/Annotation/' + self.name + '/fgs -complete=0 -train=./' + self.error_model)
        mtools.run_command(bashCommand)
    
    def annotation(self):
        diamond = DIAMOND(threads = '6',
                          db = self.db,
                          out = self.out_dir + '/Annotation/' + self.name + '/aligned.blast',
                          query = self.out_dir + '/Annotation/' + self.name + '/fgs.faa',
                          un = self.out_dir + '/Annotation/' + self.name + '/unaligned.fasta',
                          unal = '1',
                          max_target_seqs = '1')
        
        if not os.path.isfile(self.db):
            diamond.set_database(self.db, self.out_dir + '../Databases/DIAMOND')
            
        diamond.run()
    
    def uniprot_request(self, ids, original_database = 'ACC+ID', database_destination = ''):
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
            'format':'tab',
            'query':'+OR+'.join(['accession:'+acc for acc in ids]),
            'columns':'entry_name,id,ec,lineage(SUPERKINGDOM),lineage(PHYLUM),lineage(CLASS),lineage(ORDER),lineage(FAMILY),lineage(GENUS),lineage(SPECIES),pathway,protein names,database(KEGG)'
        }
        if database_destination != '':
            params['to'] = database_destination
        try:
            return http_post(BASE_URL, params=params).decode('utf-8')
        except:
            try:
                return http_post(BASE_URL, params=params).decode('utf-8')
            except:
                return ''
            
    def get_uniprot_information(self, ids, original_database = 'ACC+ID', database_destination = '', chunk = 1000):
        pbar = ProgressBar()
        result = pd.DataFrame()
        print('Retrieving UniProt information from ' + str(len(ids)) + ' IDs.')
        for i in pbar(range(0, len(ids), chunk)):
            j = i + chunk if i + chunk < len(ids) else len(ids)
            data = self.uniprot_request(ids[i:j], original_database, database_destination)
            time.sleep(60)
            if len(data) > 0:
                uniprot_info = pd.read_csv(StringIO(data), sep = '\t')
                result = pd.concat([result, uniprot_info])
        return result
    
    def recursive_uniprot_information(self, blast, output, max_iter = 5):
        if os.path.isfile(output):
            result = pd.read_csv(output, sep = '\t')
        else:
            result = pd.DataFrame()
        ids_done = set(list(result['Entry']))
        all_ids = set([ide.split('|')[1] for ide in DIAMOND(out = blast).parse_result()['sseqid']])
        i = 0
        ids_missing_output = '/'.join(output.split('/'))[:-1]
        print('Checking which IDs are missing information.')
        pbar = ProgressBar()
        ids = [ide for ide in pbar(all_ids) if ide not in ids_done]
        while len(ids_done) < len(all_ids) and i < max_iter:
            print('Information already gathered for ' + str(len(ids_done)) + 
                  ' ids. Still missing for ' + str(len(all_ids) - len(ids_done)) + '.')
            uniprotinfo = self.get_uniprot_information(ids)
            result = pd.concat([result, uniprotinfo])
            ids_done = set(list(result['Entry name']))
            print('Checking which IDs are missing information.')
            pbar = ProgressBar()
            ids = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
            i += 1
        if i < max_iter:
            print('Results for all IDs are available at ' + output)
        else:
            handler = open(ids_missing_output, 'w')
            handler.write('\n'.join(ids))
            print('Maximum iterations were made. Results related to ' + str(len(ids)) + 
                  ' IDs were not obtained. IDs with missing information are available' +
                  ' at ' + ids_missing_output + ' and information obtained is available' +
                  ' at ' + output)
        result.to_csv(output, sep = '\t')
    
    def split(self, pathway):
        pathway = pathway.split('. ')
        return [path for path in pathway if path != '']
    
    def split_ec(self, ec):
        ec = ec.split('; ')
        return [ec_number for ec_number in ec if ec_number != '']
    
    def using_repeat(self, df):
        import numpy as np
        import pandas as pd
        lens = [len(item) for item in df['Pathway']]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary["Pathway"] = np.concatenate(df['Pathway'].values)
        return pd.DataFrame(dictionary) 
    
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
            
    def info_with_coverage(self, ids, blast, uniprotinfo, assembler, mg = '', mg_reverse = ''):
        if assembler == 'metaspades':
            self.info_with_coverage_metaspades(blast, uniprotinfo)
        elif assembler == 'megahit':
            self.contigs_readcounts(output + '/Assembly/final.contigs.fa', output + '/Analysis/contigs.readcounts', 
                                    mg, mg_reverse = '')
            self.info_with_coverage_megahit(blast, output + '/Analysis/contigs.readcounts', output)
        else:
            print('Assembler not valid!')
            
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
    def run_rpsblast(self, fasta, output, cog = 'DomainIdentification/Cog'):
        mtools.run_command('rpsblast -i ' + fasta + ' -d ' + cog + ' -o ' + output + ' -m 8')
        
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
    def annotate_cogs(self, blast, output, cddid = 'DomainIdentification/cddid.tbl',
                      fun = 'DomainIdentification/fun.txt', 
                      whog = 'DomainIdentification/whog', 
                      name_of_cdd2cog = 'DomainIdentification/cdd2cog.pl'):
        mtools.run_command('perl ' + name_of_cdd2cog + ' -r ' + blast + ' -c ' + cddid + ' -f ' + fun + ' -w ' + whog)
        os.rename('results', output + '/results')
        
    '''
    Input: the output from cdd2go, a blast file with CDD and COG annotations (cogblast)
        the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt (fun)
    Output: a pandas.DataFrame with the functional categories intrisic levels 
        reorganized into corresponding columns
    '''      
    def organize_cdd_blast(self, cogblast, fun):
        cogblast = self.parse_cogblast(cogblast)
        cog_relation = self.parse_fun(fun)
        data = [cog_relation[functional_category] for functional_category in cogblast['functional categories']]
        result = pd.DataFrame(data)
        result.columns = ['functional categories', 'general functional categories']
        result = pd.concat([result[['general functional categories', 'functional categories']], 
                            cogblast[['COG protein description','qseqid']]], axis = 1)
        return result
    
    '''
    Input: the output from cdd2go, a blast file with CDD and COG annotations (cogblast)
        the fun.txt file available at ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt (fun)
    Output: an Excel with the COG identifications counted for krona plotting
    '''      
    def write_cogblast(self, cogblast, output):
        self.organize_cdd_blast(cogblast)
        cogblast = cogblast.groupby(cogblast.columns.tolist()).size().reset_index().rename(columns={0:'count'})
        cogblast.to_excel(output, index = False)
        
    '''
    Input: name of cddblast to parse
    Output: pandas.DataFrame object
    '''      
    def parse_cogblast(self, cogblast):
        cogblast = pd.read_csv(cogblast, header=None, skiprows = 1, sep = '\t', low_memory=False)
        cogblast = cogblast[list(range(0,14))+[18]]                     #several extra columns are produced because of bad formatting
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
        cog_blast: name of the COG annotation file
        fun: name of the fun.txt file
        split_pathways: boolean, if MOSCA should split the information from the
        Pathway column provided by UniProt mapping
        readcounts_matrix: name of the file containing protein expression
    Output: 
        pandas.DataFrame with all the information
    '''
    def join_reports(self, blast, uniprotinfo, cog_blast, fun, split_pathways = False,
                     readcounts_matrix = None):
        result = DIAMOND(out = blast).parse_result()
        result.index = result.qseqid
        if readcounts_matrix is not None:
            readcountsdf = pd.read_csv(readcounts_matrix, sep = '\t', header = None)[:-5]
            result = pd.merge(result, readcountsdf, left_index = True, right_index = True)
        else:
            result['Abundance'] = [ide.split('_')[5] for ide in result.qseqid]
        result.index = [ide.split('|')[1] for ide in result.sseqid]
        uniprotinfo = pd.read_csv(uniprotinfo, sep= '\t')
        uniprotinfo = uniprotinfo.drop_duplicates()
        if split_pathways == True:
            funcdf = uniprotinfo[uniprotinfo.Pathway.notnull()][['Query','Pathway']]
            funcdf.Pathway = funcdf.Pathway.apply(self.split)
            funcdf = self.using_repeat(funcdf)
            pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) 
                                        for path in funcdf.Pathway], index = funcdf.index)
            pathways.columns = ['Superpathway','Pathway','Subpathway']
            del funcdf['Pathway']; del uniprotinfo['Pathway']
            funcdf = pd.concat([pathways, funcdf], axis = 1)
            uniprotinfo = pd.merge(uniprotinfo, funcdf, on = ['Query'], how = 'outer')
        result = pd.merge(result, uniprotinfo, left_index=True, right_on = ['Entry'], how = 'outer')
        cog_blast = self.organize_cdd_blast(cog_blast, fun)
        result = pd.merge(result, cog_blast, left_on = ['qseqid'], right_on = ['qseqid'], how = 'outer')
        result = result[['Entry','Taxonomic lineage (SUPERKINGDOM)',
                        'Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)',
                        'Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)',
                        'Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)',
                        'Pathway','Protein names','EC number','Cross-reference (KEGG)',
                        'functional categories', 'general functional categories',
                        'COG protein description','Abundance']]
        result.columns = ['UniProt ID','Taxonomic lineage (SUPERKINGDOM)',
                        'Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)',
                        'Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)',
                        'Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)',
                        'Pathway','Protein name','EC number','Cross-reference (KEGG)',
                        'COG general functional category','COG functional category', 
                        'COG protein description','Abundance']
        result = result[result['ORF name'].notnull()]
        return result
        
    def run(self):
        self.gene_calling()
        self.annotation()
        blast = DIAMOND(out = self.out_dir + '/Annotation/' + self.name + '/aligned.blast').parse_result()
        ids = list(set([ide.split('|')[1] for ide in blast['sseqid'][blast.sseqid != '*']]))
        uniprot_info = self.get_uniprot_information(ids)
        uniprot_info.to_csv(self.out_dir + '/Annotation/' + self.name + '/uniprot.info', sep = '\t')
        self.recursive_uniprot_information(self.out_dir + '/Annotation/' + self.name + '/aligned.blast', 
                                           self.out_dir + '/Annotation/' + self.name + '/uniprot.info')
        self.run_rpsblast(self.out_dir + '/Annotation/' + self.name + '/fgs.faa',
                          )
        
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
    
if __name__ == '__main__':
    
    annotater = Annotater()
    
    blast = DIAMOND(out = 'MOSCAfinal/Annotation/aligned.blast').parse_result()
    
    ids = list(set([ide.split('|')[1] for ide in blast.sseqid if ide != '*']))
    
    uniprotinfo = annotater.get_uniprot_information(ids)
    
    uniprotinfo.to_csv('MOSCAfinal/Annotation/uniprot.info', sep = '\t')
    
    import shutil
    shutil.copy('MOSCAfinal/Annotation/uniprot.info', 'MOSCAfinal/Annotation/uniprot.info1')
    
    
    annotater.recursive_uniprot_information('MOSCAfinal/Annotation/aligned.blast', 'MOSCAfinal/Annotation/uniprot.info')
    '''       
    report = annotater.join_reports('MOSCAfinal/Annotation/aligned.blast',
                                       'MOSCAfinal/Annotation/uniprot.info',
                                       'MOSCAfinal/Annotation/' + sample + '/cogs.tsv',
                                       'DomainIdentification/fun.txt')
    '''
    '''
    writer = pd.ExcelWriter('MOSCAfinal/Annotation/all_info.xlsx', engine='xlsxwriter')
    
    for sample in ['4478-DNA-S1613-MiSeqKapa','4478-DNA-S1616-MiSeqKapa','4478-DNA-S1618-MiSeqKapa']:
        print(sample)
        report = annotater.join_reports('MOSCAfinal/Annotation/aligned.blast',
                                       'MOSCAfinal/Annotation/uniprot.info',
                                       'MOSCAfinal/Annotation/' + sample + '/cogs.tsv',
                                       'DomainIdentification/fun.txt')
    
        print(len(report))
        
        print(report.head())
        
        if len(report) > 1048576:
            i = 0
            k = 0
            while i + 1048576 < len(report):
                j = i + 1048576
                report.iloc[i:j].to_excel(writer, sheet_name = sample + '_' + str(k), index = False)
                i += 1048576
                k += 1
            report.iloc[i:].to_excel(writer, sheet_name = sample + '_' + str(k), index = False)
        else:
            report.to_excel(writer, sheet_name = sample, index = False)
    writer.save()
    '''
