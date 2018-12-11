# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package for retrieval of UniProt 
information and Differential Expression analysis

By João Sequeira

Sep 2017
"""

from mosca_tools import MoscaTools
import pandas as pd
import numpy as np
import glob
from annotation import Annotater

mtools = MoscaTools()

class MetaTranscriptomicsAnalyser:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
            
    def generate_mg_index(self, reference, index_prefix):
        mtools.run_command('bowtie2-build ' + reference + ' ' + index_prefix)
    
    def run_alignment(self, reads, index_prefix, sam):  
        mtools.run_command('bowtie2 -a -x ' + index_prefix + ' -1 ' + reads[0] + ' -2 ' + reads[1] + ' -S ' + sam)

    def run_htseq_count(self, sam, gff, output):
        mtools.run_command('htseq-count -i Name ' + sam + ' ' + gff, file = output)
            
    '''
    generates an expression column file
    '''
    def readcounts_file(self):
        self.generate_mg_index(self.contigs, self.out_dir + '/Analysis/' + self.mt + '_index')
        mtools.build_gff(self.out_dir + '/Annotation/' + self.mg + '/aligned.blast', self.out_dir + '/Analysis/' + self.mt + '.gff')
        self.run_alignment([self.out_dir + '/Preprocess/SortMeRNA/' + self.mt + '_' + fr + '.fastq' for fr in ['forward','reverse']],
                           self.out_dir + '/Analysis/' + self.mt + '_index', self.out_dir + '/Analysis/' + self.mt + '.sam')
        self.run_htseq_count(self.out_dir + '/Analysis/' + self.mt + '.sam', self.out_dir + '/Analysis/' + self.mt + '.gff',
                             self.out_dir + '/Analysis/' + self.mt + '.readcounts')
        
    '''
    input: 
        directory that contains the .readcounts files from htseq-count (readcounts_dir)
        names of columns of final file / names of samples (header)
    output: 
        merged expression matrix name (output)
    '''
    def merge_readcounts(self, readcounts_dir, header, output):
        files = glob.glob(readcounts_dir + '/*.readcounts')
        expression_matrix = pd.DataFrame()
        for file in files:
            df = pd.read_csv(file, sep='\t', index_col = 0, header = None)
            df.columns = [file.split('/')[-1].rstrip('.readcounts')]
            expression_matrix = pd.merge(expression_matrix, df, how = 'outer', 
                                         left_index = True, right_index = True)
        expression_matrix = expression_matrix[1:-5]                             #remove non identified proteins and the metrics at the end
        expression_matrix = expression_matrix.fillna(value=0).astype(int)       #set not identified proteins expression to 0, and convert floats, since DeSEQ2 only accepts integers
        expression_matrix.index.name = 'geneid'
        expression_matrix.columns = header
        expression_matrix.to_csv(output, sep = ' ')
    
    '''
    input: readcounts file name to change index to just COGs
        blast COG annotated file name (cog_blast)
    output: name of readcounts file to output with just COG names
    '''
    def readcounts2justcog(self, readcounts, cog_blast, output, header = None):
        pass
        
    '''
    input: readcounts file name to change index to just COGs
        blast UniProt annotated file name (blast)
        blast COG annotated file name (cog_blast)
    output: name of readcounts file to output with UniProt and COG names
    '''
    def readcounts2cogifneeded(self, readcounts, blast, cog_blast, output, header = None):
        pass
    '''
    input: readcounts file name to change index to just COGs
        blast UniProt annotated file name (blast)
        blast COG annotated file name (cog_blast)
    output: name of readcounts file to output with just Uniprot names (rest will be *)
    '''
    def readcounts2justuniprot(self, readcounts, blast, output, header = None):
        pass
    
    '''
    input: name of readcounts file concerning just ONE sample
        name of uniprotinfo file with taxonomic and functional information
    output: base name for excel files formated for krona plotting
    '''
    def readcounts2krona(self, readcounts, uniprotinfo, output):
        annotater = Annotater()
        readcountsdf = pd.read_csv(readcounts, sep = '\t', header = None, index_col = 0)[:-5]          #remove last 5 lines of htseq-count report 
        readcountsdf.columns = ['Expression']
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t', index_col = 0)
        uniprotinfo = uniprotinfo.drop_duplicates()
        joined = uniprotinfo.join(readcountsdf,how='inner')
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)','Expression']
        taxdf = joined[tax_columns].groupby(tax_columns[:-1])['Expression'].sum().reset_index()
        taxdf.to_excel(output + '_taxonomic_krona.xlsx', index = False)
        func_columns = ['Pathway','Protein names','EC number']
        funcdf = joined[joined.Pathway.notnull()][func_columns + ['Expression']]
        funcdf.Pathway = funcdf.Pathway.apply(annotater.split)
        funcdf = annotater.using_repeat(funcdf)
        pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) for path in 
                                 funcdf.Pathway], index = funcdf.index)
        pathways.columns = ['Superpathway','Pathway','Subpathway']
        del funcdf['Pathway']
        funcdf = pd.concat([pathways, funcdf], axis = 1)
        funcdf = funcdf[['Superpathway','Pathway','Subpathway','EC number',
                         'Protein names','Expression']]
        funcdf = funcdf.groupby(funcdf.columns.tolist()[:-1])['Expression'].sum().reset_index()
        funcdf.to_excel(output + '_functional_krona.xlsx', index = False)
        
    def differential_analysis(self, expression_matrix, output, conditions):
        conditions = str(conditions)[1:-1].replace("'",'"')
        print('conditions:',conditions)
        commands = ['library(DESeq2)','library(pheatmap)',
                    'total = read.table("' + output + '/all_experiments.readcounts", h=T, row.names=1)',
                    'condition <- factor(c(' + conditions + '))','total <- total[ rowSums(total) > 1, ]',
                    'cd=data.frame(c(' + conditions + '))',
                    'colnames(cd)[1]="condition"','rownames(cd)=colnames(total)',
                    'dds <- DESeqDataSetFromMatrix(countData = total,colData = cd, design = ~ condition)',
                    'dds <- DESeq(dds)',
                    'res <- results(dds)','resOrdered <- res[order(res$padj),]',
                    'jpeg("' + output + '/ma.jpeg")',
                    'plotMA(res, main="DESeq2", ylim=c(-2,2))','dev.off()',
                    'jpeg("' + output + '/counts.jpeg")', 
                    'plotCounts(dds, gene=which.min(res$padj), intgroup="condition")',
                    'dev.off()', 'write.csv(as.data.frame(resOrdered), file="' + output + '/condition_treated_results.csv")',
                    'vsd <- varianceStabilizingTransformation(dds, blind=FALSE)',
                    'select=rownames(head(resOrdered,20))','vsd.counts = assay(vsd)[select,]',
                    'df <- as.data.frame(colData(dds)[,c("condition")])',
                    'jpeg("' + output + '/gene_expression.jpeg")','pheatmap(vsd.counts)',
                    'dev.off()','sampleDists <- dist(t(assay(vsd)))',
                    'sampleDistMatrix <- as.matrix(sampleDists)',
                    'rownames(sampleDistMatrix) <- dds$condition',
                    'colnames(sampleDistMatrix) <- NULL','library(RColorBrewer)',
                    'colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)',
                    'jpeg("' + output + '/sample_distances.jpeg")',
                    'pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)',
                    'dev.off()','jpeg("' + output + '/pca.jpeg")',
                    'plotPCA(vsd, intgroup=c("condition"))','dev.off()']
        open(output + '/de_analysis.R', 'w').write('\n'.join(commands))
        mtools.run_command('Rscript ' + output + '/de_analysis.r')
        
if __name__ == '__main__':
    relation = {'4478-R1-1-MiSeqKapa':'4478-DNA-S1613-MiSeqKapa','4478-R2-1-MiSeqKapa':'4478-DNA-S1613-MiSeqKapa',
                '4478-R3-1-MiSeqKapa':'4478-DNA-S1616-MiSeqKapa','4478-R4-1-MiSeqKapa':'4478-DNA-S1618-MiSeqKapa'}
    
    for mt in ['4478-R1-1-MiSeqKapa', '4478-R3-1-MiSeqKapa']:
        
        toolbox = MetaTranscriptomicsAnalyser(assembler = 'metaspades')
        
        toolbox.readcounts2krona('MOSCAfinal/Analysis/' + mt + '.readcounts',
                                 'MOSCAfinal/Annotation/' + relation[mt] + '/uniprot.info',
                                 'MOSCAfinal/Analysis/' + mt)