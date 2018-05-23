# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:02:36 2018

@author: Asus
"""

class Archive
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
        intvalues = gis.values.astype(int)
        ngis = pd.DataFrame(intvalues)
        taxons = self.lineages(taxids)
        taxons = taxons[taxons.index.notnull()]
        systems = self.pathways(ngis)
        self.write_df(taxons, out + '/taxonomy.xlsx')
        self.write_df(systems, out + '/systems.xlsx')
        
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
        
    def build_gtf(self, blast, db, output):
        import re, csv
        pattern = re.compile('(.+cov_.+?)_(.+?)_(.+?)_(\+|-)')
        parts = [re.search(pattern, q).groups() for q in blast.qseqid]
        dimension = len(parts)
        lines = [[line[0] for line in parts]]
        lines.append([db for i in range(dimension)])
        lines.append(['exon' for i in range(dimension)])
        lines.append([line[1] for line in parts])
        lines.append([line[2] for line in parts])
        lines.append([score for score in blast.evalue])
        lines.append(lines.append([line[3] for line in parts]))
        lines.append(['.' for i in range(dimension)])
        lines.append(['gene_name "' + name.split('|') if name != '*' else name + '"' for name in blast.sseqid])
        gtf = pd.DataFrame(lines)
        gtf = gtf.transpose()
        gtf.write_csv(output + '/gtf.gtf', sep = '\t', header = False, quoting=csv.QUOTE_NONE, index=False)
        

        
    def generate_sam(self, contigs, output, reads, paired):
        first_command = 'bowtie2-build ' + contigs + ' ' + output + '/idx'
        second_command = {'PE':'bowtie2 -x ' + output + '/idx -1 ' + reads[0] + ' -2 ' + reads[1] + ' -S ' + output + '/sam.sam',
                          'SE':'bowtie2 -x ' + output + '/idx -U ' + reads[0] + ' -S ' + output + '/sam.sam'}
        commands = (first_command, second_command[paired])
        for command in commands:
            self.run(command)
            
    def readcounts(self, sam, gtf, output):
        command = 'htseq-count -i gene_name ' + sam + ' ' + gtf + ' > ' + output + '/readcounts.readcounts'
        self.run(command)
        
    def r_analysis(self, readcounts, output):
        cat_command = 'echo "genename'
        for file in readcounts: cat_command += ' ' + file
        cat_command += '" > ' + output + '/expression.tab'
        paste_command = 'paste ' + self.working_dir + '/Analysis/*.readcounts | cut -f1,2,4,6,8,10,12 >> ' + output + '/expression.tab
    
    def differential_expression(self, assembler, reads):
        contigs = {'metaspades':'contigs.fasta', 'megahit':'final.contigs.fa'}
        self.generate_sam(contigs[assembler], self.working_dir + '/Analysis', reads, self.paired)
        script = self.de_script()
        self.build_gtf(self.working_dir + '/Annotation/aligned.blast', self.db, self.working_dir + '/Analysis')
        self.readcounts(self.working_dir + '/Analysis/sam.sam', self.working_dir + '/Analysis/gtf.gtf', self.working_dir + '/Analysis)
        self.r_script()
        
    def run(self, bashCommand):
        print(bashCommand)
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        return output, error
    
    def parse_uniprot_tab(self, file):
        df = pd.DataFrame.from_csv(file, sep = '\t')
        tax_columns = ['Taxonomic lineage (' + tax + ')' for tax in ['SUPERKINGDOM', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']]