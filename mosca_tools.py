# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:17:58 2017

@author: Asus
"""
'''
# ssh -X bridge_root@193.137.11.212
sudo apt-get install screen
'''
import subprocess
from progressbar import ProgressBar
from io import StringIO
import pandas as pd
import glob
import re

class MoscaTools:
    
    def run_command(self, bashCommand, file = '', mode = 'w'):
        print(bashCommand)
        if file == '':
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        else:
            handler = open(file, mode)
            process = subprocess.Popen(bashCommand.split(), stdout=handler)
        output, error = process.communicate()
    
    def build_gff(self, aligned, output):
        import pandas as pd
        gff = pd.DataFrame()
        diamond = DIAMOND(out = aligned).parse_result()
        parts = [qid.split('_') for qid in diamond.qseqid]
        preid = [part[1] for part in parts]
        node = 1
        j = 1
        ids = []
        for i in preid:
            if i == node:
                ids.append('seq' + str(i) + '_' + str(j))
                j += 1
            else:
                node = i
                ids.append('seq' + str(i) + '_' + str(j))
                j = 1
        gff["seqid"] = ['NODE_' + part[1] + '_' + part[2] + '_' + part[3] + '_' + part[4] + '_' + part[5] for part in parts]
        size = gff.size
        gff["source"] = ['UniProtKB' for i in range(size)]
        gff["type"] = ['exon' for i in range(size)]
        gff["start"] = [part[6] for part in parts]
        gff["end"] = [part[7] for part in parts]
        gff["score"] = diamond.evalue
        gff["strand"] = [part[8] for part in parts]
        gff["phase"] = ['.' for i in range(size)]
        gff["ID"] = [ide.split('|')[2] if ide != '*' else ide for ide in diamond.sseqid]
        gff["Name"] = diamond.qseqid
        gff["attributes"] = ['gene_id=' + i.Name + ';Name=' + i.ID for i in gff.itertuples()]
        del gff['ID']; del gff['Name']
        gff.to_csv(output + '.gff', sep = '\t', index=False, header=False)
        return gff
    
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
        import subprocess
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
        self.run_command('bash MOSCA/merge-paired-reads.sh ' + file1 + ' ' + file2 + ' ' + output)
    
    def divide_fq(self, file, output1, output2):
        self.run_command('bash MOSCA/unmerge-paired-reads.sh ' + file + '.fastq ' + output1 + ' ' + output2)
        
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
            
    def set_to_uniprotID(self, fasta, aligned, output):
        handler = DIAMOND(out = aligned)
        result = handler.parse_result(aligned)
        sequences = self.parse_fasta(fasta)
        keys = list(sequences.keys())
        final = list()
        for name in keys:
            try:
                newname = str(result[result.qseqid == name]['sseqid'].item())
                final.append((newname, sequences[name]))
            except:
                print(result[result.qseqid == name]['sseqid'])
        with open(output,'w') as f:
            for double in final:
                f.write('>' + double[0] + '\n' + double[1] + '\n')
                
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
    resultdf.loc[resultdf['lineage(FAMILY)'].str.contains('Chloroflexaceae'),'lineage(FAMILY)'] = 'Chloroflexaceae'
    resultdf.loc[resultdf['lineage(SPECIES)'].str.contains('Methanosarcina mazei'),'lineage(SPECIES)'] = 'Methanosarcina mazei'
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
    resultdf.columns = ['id','ec','lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)','lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','pathway','protein names','coverage']
    resultdf = resultdf[['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)','lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','coverage']]
    resultdf = resultdf.groupby(resultdf.columns.tolist()[:-1])['coverage'].sum().reset_index()
    realitydf = pd.read_csv(reality, sep = '\t', header = None)
    realitydf.columns = ['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                'lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','abundance']
    realitydf['lineage(SPECIES)'] = realitydf['lineage(GENUS)'] + ' ' + realitydf['lineage(SPECIES)']
    total = float(resultdf['coverage'].sum())
    print('Total is ' + str(total))
    resultdf['coverage'] /= (total/100)      #set as percentage
    totalreal = float(realitydf['abundance'].sum())
    realitydf['abundance'] /= (totalreal/100)
    #ammendments on the result because of some uniprot facts
    resultdf.loc[resultdf['lineage(FAMILY)'].str.contains('Chloroflexaceae'),'lineage(FAMILY)'] = 'Chloroflexaceae'
    resultdf.loc[resultdf['lineage(SPECIES)'].str.contains('Methanosarcina mazei'),'lineage(SPECIES)'] = 'Methanosarcina mazei'
    resultdf.loc[resultdf['lineage(GENUS)'].str.contains('Methanothrix'),'lineage(GENUS)'] = 'Methanosaeta'
    for level in ['lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                'lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)']:
        print('Level: ' + level)
        factors = list(realitydf[level].unique())
        value = 0
        others = 100
        for factor in factors:
            print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['coverage'].sum())) + '%')
            others -= float(resultdf[resultdf[level] == factor]['coverage'].sum())
            value += (float(resultdf[resultdf[level] == factor]['coverage'].sum()) - abs(float(resultdf[resultdf[level] == factor]['coverage'].sum()) - realitydf.loc[realitydf[level] == factor]['abundance'].sum()))
        if level == 'lineage(SUPERKINGDOM)':
            for factor in ['Eukaryota','Viruses']:
                print('Percentage of ' + factor + ': ' + str(float(resultdf[resultdf[level] == factor]['coverage'].sum())) + '%')
                others -= float(resultdf[resultdf[level] == factor]['coverage'].sum())
        print('Others:' + str(others) + '%')
        print('Correctness: ' + str(value) + '%')
        
for value in simulateddf1.index.get_values():
    print(value)
    found = False
    partdf = metaspadestaxdf1[metaspadestaxdf1.index.str.contains(value)]
    print(partdf.head())
    for i in range(len(metaspadestaxdf1)):
        print(metaspadestaxdf1[~(metaspadestaxdf1==simulateddf1).all(axis=1)])
        
    df[~(df==df1).all(axis=1)]    
    
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
    lens = [len(item) for item in df['pathway']]
    dictionary = dict()
    for column in df.columns:
        dictionary[column] = np.repeat(df[column].values,lens)
    dictionary["pathway"] = np.concatenate(df['pathway'].values)
    return pd.DataFrame(dictionary)

def partition_pathway(pathway):
    return [(path[0].split('; ') + [np.nan] * (3 - len(path[0].split('; ')))) for path in pathway]
            
def organize_pathway(result, name):
    import pandas as pd    
    import numpy as np
    resultdf = pd.read_csv(result, sep = '\t', header = None)
    resultdf.columns = ['id','ec','lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)',
                        'lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','pathway','protein names']
    resultdf = resultdf[resultdf.pathway.notnull()]
    resultdf.pathway = resultdf.pathway.apply(split)
    resultdf = using_repeat(resultdf)
    pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) for path in 
                             resultdf.pathway], index = resultdf.index)
    resultdf = resultdf[['protein names','id','ec','lineage(GENUS)']]
    pathways.columns = ['superpathway','pathway','subpathway']
    resultdf = pd.concat([pathways, resultdf], axis = 1)
    resultdf = resultdf.groupby(resultdf.columns.tolist()).size().reset_index().rename(columns={0:name})
    return resultdf.set_index(resultdf.columns.tolist()[:-1])

def make_contigs_gff(file, output):
    pbar = ProgressBar()
    print('Building gff file from Megahit contigs at',file,'to',output)
    mt = MoscaTools()
    contigs = mt.parse_fasta(file)
    gff = list()
    for name in pbar(contigs.keys()):
        gff.append([name, 'Megahit', 'exon', '1', str(len(contigs[name])), '.', '.', '.', 'Name=' + name])
    print(gff)
    gff = pd.DataFrame(gff)
    gff.to_csv(output + '/gff.gff', sep = '\t', index=False, header=False)
    
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
    