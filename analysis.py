# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package for retrieval of UniProt 
information and Differential Expression analysis

By João Sequeira

Sep 2017
"""

from mosca_tools import MoscaTools
from progressbar import ProgressBar
from diamond import DIAMOND
import pandas as pd

mtools = MoscaTools()

class MetaTranscriptomicsAnalyser:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
            
    def build_gff(self, blast, output):
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
    def run(self):
        self.build_gff(self.out_dir + '/Annotation/aligned.blast', self.out_dir + '/Analysis')
        reference = {'metaspades':'contigs.fasta','megahit':'final.contigs.fa'}
        self.run_alignment(self.mt, self.out_dir + '/Assembly/' + reference[self.assembler], self.out_dir + '/Analysis')
        self.run_R(self.out_dir + '/Analysis', self.out_dir + '/Analysis', self.conditions)