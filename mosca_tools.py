# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:17:58 2017

@author: Asus
"""

import subprocess

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
        lines = [line.rstrip('\n') for line in open(file)]
        i = 0
        sequences = dict()
        print(len(lines))
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
                
    def uniprot_mapping(self, ids, max_retry = 3):
        import requests, time
        
        BASE_URL = 'http://www.uniprot.org/uniprot/'
        
        params = {
            'from':'ACC+ID',
            'to':'ACC+ID',
            'format':'tab',
            'query':'+OR+'.join(['accession:'+acc for acc in ids]),
            'columns':'id,ec,lineage(SUPERKINGDOM),lineage(PHYLUM),lineage(CLASS),lineage(ORDER),lineage(FAMILY),lineage(GENUS),lineage(SPECIES),pathway,protein names'
            }
        response = requests.post(BASE_URL, params=params)
        cnt=0
        max_retry = max_retry
        while cnt < max_retry and response.status_code != requests.codes.ok:
            cnt += 1
            time.sleep(60)
            response = requests.post(BASE_URL, params=params)
        return [match.split('\t') for match in response.content.decode('utf-8').split('\n')[1:-1]]
        
    def chunky(self, ids, step = 100, max_retry = 3):
        result = list()
        for i in range(0, len(ids), step):
            query = ids[i: i + step]
            maped = self.uniprot_mapping(query, max_retry)
            result += self.uniprot_mapping(query, max_retry)
        return pd.DataFrame([row for row in result], index = [row[0] for row in result], columns = ['id','ec','lineage(SUPERKINGDOM)','lineage(PHYLUM)','lineage(CLASS)','lineage(ORDER)','lineage(FAMILY)','lineage(GENUS)','lineage(SPECIES)','pathway','protein names'])
    