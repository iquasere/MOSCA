#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 00:23:51 2017

@author: sequeira
"""

import pandas as pd
import numpy as np
import subprocess

class Analysing:
    def __init__(self, **kwargs):
        self.__dict__ = kwargs
        if hasattr(self, 'accession2taxid'):
            self.hdf = self.accession2taxid + '.h5'
    
    def iter_read_csv(self, chunksize, **kwargs):
        for chunk in pd.read_csv(self.accession2taxid, chunksize=chunksize, low_memory=True, **kwargs):
            yield chunk
    
    def write_to_hdf_store(self, columns=None):
        hdf = pd.HDFStore(self.hdf, complevel=3)
        generator = self.iter_read_csv(self.chunksize, delimiter="\t")
        first = next(generator)
        original_cols = first.columns
        if columns is not None:
            first.columns = columns
        hdf.put('all', first, format='table', append=True, chunksize=self.chunksize, data_columns=original_cols if columns is None else columns, index=False)
        for chunk in generator:
            if columns is not None:
                chunk.columns = columns
            hdf.append('all', chunk, format='table', append=True, chunksize=self.chunksize, data_columns=original_cols if columns is None else columns, index=False)
        if columns is not None:
            hdf.create_table_index('all', columns=columns)
        else:
            hdf.create_table_index('all', columns=original_cols)
            
    def ids(self, accession):
        handler = pd.HDFStore(self.hdf, key='all', chunksize = self.chunksize)
        return handler.select('all', 'accession_version=' + str(list(accession)), chunksize = self.chunksize)
    
    def lineages(self, taxids):
        df = pd.read_csv(self.lineages)
        df.index = df.tax_id
        return pd.concat([taxids,df.loc[taxids.index,['superkingdom','phylum','class','order','family','genus','species']]], axis = 1)
    
    def pathways(self, gis):
        df = pd.read_csv(self.systems, encoding = "ISO-8859-1", names = ['bsid','source','accession','name','type','taxonomic_scope','taxid','description'], sep = '\t')
        df.index = df.gi
        return pd.concat([gis,df.loc[gis, ['name','type']]])

    def krona_plots(self, accession, out):
        if self.contigs == True:
            accession = accession[accession.values != '*']
            prelation = self.ids(accession)
        else:
            accession = accession.groupby(accession.sseqid).size().reset_index(name='count')
            prelation = self.ids(accession.sseqid)
        print(len(accession))
        for i in prelation:
            relation = i
        relation.index = relation.accession_version
        #unique, counts = np.unique(relation, return_counts=True)
        taxids, gis = relation.loc[accession,'taxid'], relation.loc[accession,'gi']
        taxids, gis = taxids.groupby(taxids.values).size().reset_index(name='counting'), gis.groupby(gis.values).size().reset_index(name='counting')
        taxids, gis = taxids[taxids.counting > (0.001 * sum(taxids.counting))], gis[gis.counting > (0.00005 * sum(gis.counting))]
        #intvalues = gis.values.astype(int)
        #ngis = pd.DataFrame(intvalues)
        taxons = self.lineages(taxids)
        taxons = taxons[taxons.index.notnull()]
        #systems = self.pathways(ngis)
        self.write_df(taxons, out + '/taxonomy.xlsx')
        #self.write_df(systems, out + '/systems.xlsx')
        
        print("'tis done!")
        
    def write_df(self, df, out):
        print('creating excel')
        writer = pd.ExcelWriter(out, engine='xlsxwriter')
        df.to_excel(writer,'Sheet1')
        writer.save()
        print('finished csv')
        
    def test_h5_tax(self):
        self.write_to_hdf_store(['accession', 'accession_version', 'taxid', 'gi'])
        self.krona_plots(self.ids)
        
    def index_contigs(self):
        contigs = {'megahit':'/final.contigs.fa', 'metaspades':'/contigs.fasta'}
        bashCommand = 'bowtie2-build ' + self.working_dir + '/Annotation' + contigs[self.assembler] + self.working_dir + '/Assembly/bowtie2/idx'
        self.run(bashCommand)
        
    def de_matrix(self, blasts, out):
        from diamond import DIAMOND
        result = pd.DataFrame()
        for blast in blasts:
            data = DIAMOND(out = blast).parse_result()
            df = pd.DataFrame()
            df[blast.split('/')[0]] = data[data['length'] > 0].qseqid.apply(lambda x: x.split('_')[5]).astype(float)
            print(df)
            #convert to RPK
            df[blast.split('/')[0]] = df[blast.split('/')[0]] / (data[data['length'] > 0]['length'] / 1000)
            print(df)
            #convert to TPM            
            per_M = sum(df[blast.split('/')[0]]) / 1000000
            print(per_M)
            df[blast.split('/')[0]] = df[blast.split('/')[0]] / per_M
            print(df)
            df.index = data[data['length'] > 0].sseqid.apply(lambda x: x.split('|')[-1])            
            sdf = df.groupby(df.index).sum()
            result = pd.concat([result,sdf], axis = 1)
        result.fillna(0, inplace=True)
        result = result.fillna(0)
        result.to_csv(out, sep = '\t')
        
    def de_script(self, blast):
        pass
    
    def differential_expression(self):
        #http://software.broadinstitute.org/software/igv/download#binary
        #https://github.com/deweylab/RSEM
        #sudo apt-get install perl-doc for rsem-plot-transcript-wiggles
        contigs = {'metaspades':'contigs.fasta', 'megahit':'final.contigs.fa'}
        if self.paired == 'PE':
            commands = (['bwa index ' + self.project + '/Assembly/MetaSPAdes/' + contigs.fasta + '-p Assembly/idx'] +
            ['bwa aln ' + self.project + '/Assembly/idx ' + self.project + '/Preprocess/SortMeRNA/' + fr + '_rejected > Assembly/' + fr + '_rejected.sai' for fr in ['forward','reverse']] +
            ['bwa sampe ' + self.project + '/Assembly/idx ' + self.project + '/Assembly/forward_rejected.sai ' + self.project + '/Assembly/reverse_rejected.sai ' + self.project + '/Preprocess/SortMeRNA/forward_rejected ' + self.project + '/Preprocess/SortMeRNA/reverse_rejected | sam view -> ' + self.project + '.sam'])
        script = self.de_script()
        
    def run(self, bashCommand):
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error