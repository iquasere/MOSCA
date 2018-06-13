# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package for retrieval of UniProt 
information and Differential Expression analysis

By João Sequeira

Sep 2017
"""

import pandas as pd
import glob
from mosca_tools import MoscaTools
from mosca_tools import get_ids_from_ncbi_gff
from progressbar import ProgressBar
from io import StringIO
from diamond import DIAMOND

mtools = MoscaTools()

class Analyser:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs

        
    def uniprot_request(self, ids, original_database = 'ACC+ID', database_destination = ''):
        import requests
         
        BASE_URL = 'http://www.uniprot.org/uniprot/'
         
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
        
        try:
            params = {
                'from':original_database,
                'format':'tab',
                'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                'columns':'id,ec,lineage(SUPERKINGDOM),lineage(PHYLUM),lineage(CLASS),lineage(ORDER),lineage(FAMILY),lineage(GENUS),lineage(SPECIES),pathway,protein names,database(KEGG)'
            }
            if database_destination != '':
                params['to'] = database_destination
            response = http_post(BASE_URL, params=params)
            time.sleep(3)
            return response.decode('utf-8')
        except:
            print(ids[0],'to',ids[-1],'tab failed')
            time.sleep(120)
            try:
                params = {
                    'from':original_database,
                    'format':'tab',
                    'query':'+OR+'.join(['accession:'+acc for acc in ids]),
                    'columns':'id,ec,lineage(SUPERKINGDOM),lineage(PHYLUM),lineage(CLASS),lineage(ORDER),lineage(FAMILY),lineage(GENUS),lineage(SPECIES),pathway,protein names,database(KEGG)'
                    }
                if database_destination != '':
                    params['to'] = database_destination
                response = http_post(BASE_URL, params=params)
                time.sleep(3)
                return response.decode('utf-8')
            except:
                print(ids[0] + ' to ' + ids[-1] + ' tab failed')
            
    def get_uniprot_information(self, ids, original_database = 'ACC+ID', database_destination = '', chunk = 100):
        pbar = ProgressBar()
        result = pd.DataFrame()
        #print('Writing Uniprot information to ' + output)
        for i in pbar(range(0, len(ids), chunk)):
            if i + chunk < len(ids):
                data = StringIO(self.uniprot_request(ids[i:i+chunk], original_database, database_destination))
                uniprot_info = pd.read_csv(data, sep = '\t')
                result = pd.concat([result, uniprot_info])
            else:
                data = StringIO(self.uniprot_request(ids[i:len(ids)], original_database, database_destination))
                uniprot_info = pd.read_csv(data, sep = '\t')
                return pd.concat([result, uniprot_info])
            
    def build_gff(self, blast, output):
        from diamond import DIAMOND
        gff = pd.DataFrame()
        diamond = DIAMOND(out = blast).parse_result()
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
        if self.assembler == 'metaspades':
            gff["seqid"] = ['_'.join(part[:-3]) for part in parts]
        elif self.assembler == 'megahit':    #must first run this command: sed -i 's/\s.*$//' contigs.fa
            constant = parts[0][0]
            gff["seqid"] = [constant + '_' + part[1] for part in parts[:-3]]
        size = gff.size
        gff["source"] = ['UniProtKB' for i in range(size)]
        gff["type"] = ['exon' for i in range(size)]
        gff["start"] = [part[-3] for part in parts]
        gff["end"] = [part[-2] for part in parts]
        gff["score"] = diamond.evalue
        gff["strand"] = [part[-1] for part in parts]
        gff["phase"] = ['.' for i in range(size)]
        gff["ID"] = [ide.split('|')[2] if ide != '*' else ide for ide in diamond.sseqid]
        gff["Name"] = diamond.qseqid
        gff["attributes"] = ['gene_id=' + i.Name + ';Name=' + i.ID for i in gff.itertuples()]
        del gff['ID']; del gff['Name']
        gff.to_csv(output + '/gff.gff', sep = '\t', index=False, header=False)
    
    #by experiment
    def run_bowtie2(self, forward_reads, reverse_reads, reference, output, number = '', fasta = ''):   #fasta = '' or ' -f' 
        mtools.run_command('bowtie2 -a' + fasta + ' -x ' + output + '/idx -1 ' + forward_reads + ' -2 ' + reverse_reads + ' -S ' + output + '/mt' + str(number) + '.sam')

    #for all experiments
    def run_htseq_count(self, sam, gff, output, number = ''):
        mtools.run_command('htseq-count -i Name ' + sam + ' ' + gff, file = output + '/readcounts' + str(number) + '.readcounts')
        
    #generate index from contigs, perform alignment and measure expression
    def run_alignment(self, files, reference, output):   #files = [[forward_reads, reverse_reads], ...]
        number = 1
        mtools.run_command('bowtie2-build ' + reference + ' ' + output + '/idx')
        for group in files:
            print('Aligning ' + group[0] + ' and ' + group[1] + ' to ' + reference + ' as number ' + str(number))
            self.run_bowtie2(group[0], group[1], reference, output, number = number, fasta = ' -f')
            self.run_htseq_count(output + '/mt' + str(number) + '.sam', output + '/gff.gff', output, number = number)
            number += 1
        
    def run_R(self, readcounts_dir, output, conditions):
        #join readcounts files
        files = [readcounts_dir + '/readcounts' + str(i) + '.readcounts' for i in range(1,len(conditions) + 1)]
        result = pd.read_csv(files[0], header = None, sep = '\t', index_col = 0)
        for file in files[1:]:
            result = pd.concat([result, pd.read_csv(file, header = None, sep = '\t', index_col = 0)], axis = 1)
        del result.index.name
        result.columns = conditions
        result.to_csv(output + '/total.tab', sep = '\t')
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
        mtools.run_command('Rscript ' + output + '/de_analysis.R')
            
    #for multi sample
    def differential_analysis(self):
        self.build_gff(self.out_dir + '/Annotation/aligned.blast', self.out_dir + '/Analysis')
        reference = {'metaspades':'contigs.fasta','megahit':'final.contigs.fa'}
        self.run_alignment(self.mt, self.out_dir + '/Assembly/' + reference[self.assembler], self.out_dir + '/Analysis')
        self.run_R(self.out_dir + '/Analysis', self.out_dir + '/Analysis', self.conditions)

    def full_analysis(self):
        self.info_with_coverage(self.out_dir + '/Annotation/aligned.blast', self.out_dir)
        
        if self.de == True:
            self.differential_analysis()
            
    def info_with_coverage(self, blast, output):
        if self.assembler == 'metaspades':
            self.info_with_coverage_metaspades(blast, output + '/Analysis/uniprot_information_coverage.tab')
        elif self.assembler == 'megahit':
            self.contigs_readcounts(self.mt[0], self.mt[1], output + '/Assembly/final.contigs.fa', 
                                    output + '/Analysis/contigs.readcounts')
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
        
    def contigs_readcounts(self, read1, read2, contigs, output):
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
        
    def split(self, pathway):
        pathway = pathway.split('. ')
        return [path for path in pathway if path != '']
    
    def using_repeat(self, df):
        import numpy as np
        import pandas as pd
        lens = [len(item) for item in df['Pathway']]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary["Pathway"] = np.concatenate(df['Pathway'].values)
        return pd.DataFrame(dictionary)
    
    def organize_pathway(self, fullinfo):
        import pandas as pd    
        import numpy as np
        resultdf = pd.read_csv(fullinfo, sep = '\t')
        print(resultdf)
        resultdf = resultdf[resultdf.Pathway.notnull()]
        resultdf.pathway = resultdf.Pathway.apply(self.split)
        resultdf = self.using_repeat(resultdf)
        pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) for path in 
                                 resultdf.Pathway], index = resultdf.index)
        resultdf = resultdf[['Protein names','Entry','EC number','Taxonomic lineage (GENUS)']]
        pathways.columns = ['Superpathway','Pathway','Subpathway']
        del resultdf.Pathway
        resultdf = pd.concat([pathways, resultdf], axis = 1)
        print(resultdf.head())
        columns = resultdf.columns.tolist()
        columns.remove('coverage')
        resultdf = resultdf.groupby(columns)['coverage'].sum().reset_index()
        return resultdf.set_index(resultdf.columns.tolist()[:-1])

    def organize_func(self, resultdf):
        resultdf = resultdf[resultdf.Pathway.notnull()]
        resultdf.Pathway = resultdf.Pathway.str.split('\. ')
        rows = list()
        for i in range(len(resultdf)):
            for value in resultdf.iloc[i].Pathway:
                if value != '':
                    rows.append([value])
        for i in range(len(rows)):
            rows[i] = rows[i][0].split('; ')
        rows = [row + [np.nan] * (3 - len(row)) for row in rows]
        resultdf = pd.DataFrame(rows)
        resultdf.columns = ['Superpathway','Pathway','Subpathway']
        resultdf = resultdf.groupby(resultdf.columns.tolist()[:-1]).size().reset_index().rename(columns={0:'abundance'})
        return resultdf
