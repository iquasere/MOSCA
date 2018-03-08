#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 00:23:51 2017

@author: sequeira
"""

import pandas as pd
import glob
from mosca_tools import MoscaTools
mt = MoscaTools()

class Analysing:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

        
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

    def get_uniprot_information(self, ids, output, chunk):
        i = 0
        j = 0
        while j < len(ids):
            if i + chunk > len(ids):
                j = len(ids)
            else:
                j = i + chunk
            self.uniprot_request(ids[i:j], output)
            i = i + chunk
            
    def build_gff(self, blast, output):
        from diamond import DIAMOND
        gff = pd.DataFrame()
        diamond = DIAMOND(out = blast).parse_result()
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
        gff.to_csv(output + '/gff.gff', sep = '\t', index=False, header=False)
    
    #by experiment
    def run_bowtie2(self, forward_reads, reverse_reads, reference, output, number = '', fasta = ''):   #fasta = '' or ' -f' 
        mt.run_command('bowtie2 -a' + fasta + ' -x ' + output + '/idx -1 ' + forward_reads + ' -2 ' + reverse_reads + ' -S ' + output + '/mt' + str(number) + '.sam')

    #for all experiments
    def run_htseq_count(self, sam, gff, output, number = ''):
        mt.run_command('htseq-count -i Name ' + sam + ' ' + gff, file = output + '/readcounts' + str(number) + '.readcounts')
        
    #perform alignment and measure expression
    def run_alignment(self, files, reference, output):   #files = [[forward_reads, reverse_reads], ...]
        number = 1
        mt.run_command('bowtie2-build ' + reference + ' ' + output + '/idx')
        for group in files:
            print('Aligning ' + group[0] + ' and ' + group[1] + ' to ' + reference + ' as number ' + str(number))
            #self.run_bowtie2(group[0], group[1], reference, output, number = number, fasta = ' -f')
            self.run_htseq_count(output + '/mt' + str(number) + '.sam', output + '/gff.gff', output, number = number)
            number += 1
        
    def run_R(self, readcounts_dir, output, conditions):
        import os
        #join readcounts files
        if os.path.isfile(output + '/temp.readcounts'):
            os.remove(output + '/temp.readcounts')
        if os.path.isfile(output + '/total.tab'):
            os.remove(output + '/total.tab') 
        readcounts_files = glob.glob(readcounts_dir + '/*.readcounts')
        mt.run_command('paste ' + ' '.join(readcounts_files),  output + '/temp.readcounts', 'w')
        mt.run_command('cut ' + output + '/temp.readcounts -f1,2,4,6', output + '/total.tab', 'w')
        df = pd.read_csv(output + '/total.tab', sep = '\s', index_col = 0, header = None)
        df.columns = conditions
        df = df[:-5]
        df.to_csv(output + '/total.tab', sep = '\t', header = True, index = True)
        if os.path.isfile(output + '/temp.readcounts'):
            os.remove(output + '/temp.readcounts')
        #differential analysis
        conditions = str(conditions)[1:-1].replace("'",'"')
        print('conditions:',conditions)
        commands = ['library(DESeq2)','library(pheatmap)','total = read.table("' + output + '/total.tab", h=T, row.names=1)',
                    'condition <- factor(c(' + conditions + '))','total <- total[ rowSums(total) > 1, ]','cd=data.frame(c(' + conditions + '))',
                    'colnames(cd)[1]="condition"','rownames(cd)=colnames(total)',
                    'dds <- DESeqDataSetFromMatrix(countData = total,colData = cd, design = ~ condition)','dds <- DESeq(dds)',
                    'res <- results(dds)','resOrdered <- res[order(res$padj),]','jpeg("' + output + '/ma.jpeg")',
                    'plotMA(res, main="DESeq2", ylim=c(-2,2))','dev.off()','jpeg("' + output + '/counts.jpeg")', 
                    'plotCounts(dds, gene=which.min(res$padj), intgroup="condition")','dev.off()',
                    'write.csv(as.data.frame(resOrdered), file="' + output + '/condition_treated_results.csv")',
                    'vsd <- varianceStabilizingTransformation(dds, blind=FALSE)','select=rownames(head(resOrdered,20))',
                    'vsd.counts = assay(vsd)[select,]','df <- as.data.frame(colData(dds)[,c("condition")])',
                    'jpeg("' + output + '/gene_expression.jpeg")','pheatmap(vsd.counts)','dev.off()','sampleDists <- dist(t(assay(vsd)))',
                    'sampleDistMatrix <- as.matrix(sampleDists)','rownames(sampleDistMatrix) <- dds$condition',
                    'colnames(sampleDistMatrix) <- NULL','library(RColorBrewer)','colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)',
                    'jpeg("' + output + '/sample_distances.jpeg")',
                    'pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)',
                    'dev.off()','jpeg("' + output + '/pca.jpeg")','plotPCA(vsd, intgroup=c("condition"))','dev.off()']
        open(output + '/de_analysis.R', 'w').write('\n'.join(commands))
        mt.run_command('Rscript ' + output + '/de_analysis.R')
            
    #for multi sample
    def differential_analysis(self):
        self.build_gff(self.out_dir + '/Annotation/aligned.blast', self.out_dir + '/Analysis')
        self.run_alignment(self.mt, self.out_dir + '/Assembly/contigs.fasta', self.out_dir + '/Analysis')
        self.run_R(self.out_dir + '/Analysis', self.out_dir + '/Analysis', ['MT1','MT2','MT3'])

    def full_analysis(self):
        self.get_uniprot_information(self.ids, output = self.out_dir, chunk = 100)
        if self.de == True:
            self.differential_analysis()