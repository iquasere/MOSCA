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
from annotation import Annotater

mtools = MoscaTools()

class MetaTranscriptomicsAnalyser:
    
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
            
    '''
    generates an expression column file
    '''
    def readcounts_file(self):
        mtools.perform_alignment(self.contigs, self.reads, self.out_dir + '/' + self.mt, 
                                 blast = self.blast, threads = self.threads)

    '''
    input: 
        readcount_files: files from htseq-count with protein expression for each sample
        header: names of columns of final file / names of samples
        output: directory to output results
    output: 
        merged expression matrix name (output)
    '''
    def merge_readcounts(self, readcount_files, header, output):
        expression_matrix = pd.DataFrame()
        for file in readcount_files:
            df = pd.read_csv(file, sep='\t', index_col = 0, header = None)
            df.columns = [file.split('/')[-1].rstrip('.readcounts')]
            expression_matrix = pd.merge(expression_matrix, df, how = 'outer', 
                                         left_index = True, right_index = True)
        expression_matrix = expression_matrix[1:-5]                             #remove non identified proteins (*) and the metrics at the end
        expression_matrix = expression_matrix.fillna(value=0).astype(int)       #set not identified proteins expression to 0, and convert floats, since DeSEQ2 only accepts integers
        expression_matrix.index.name = 'geneid'
        expression_matrix.columns = header
        expression_matrix.to_csv(output, sep = '\t')
    
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
        
    def differential_analysis(self, readcounts, conditions, output = None):
        print('Performing differential expression analysis.')
        print('Readcounts file: ' + readcounts)
        print('Conditions: ' + ','.join(conditions))
        print('Results will be exported to: ' + output)
        mtools.run_command('Rscript MOSCA/de_analysis.R --readcounts ' + readcounts +
                         ' --conditions "' + ' '.join(conditions) + '" --output ' + output)
