"""
MOSCA's Annotation package for Gene Calling and 
Alignment of identified ORFs to UniProt database

By João Sequeira

Jun 2017
"""

from diamond import DIAMOND
from mosca_tools import MoscaTools
from uniprot_mapping import UniprotMapping
from progressbar import ProgressBar
from tqdm import tqdm
import pandas as pd
import numpy as np
import os, glob

mtools = MoscaTools()
upmap = UniprotMapping()

class Annotater:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    '''
    Input:
        file: name of input file to perform gene calling on
        output: basename of output files
        assembled: True if input is contigs, False if it are reads
        error_model: quality model to consider when input are reads
    Output:
        FragGeneScan output files will be produced with basename 'output + /fgs'
        If input is FASTQ reads (if assembled == False) a FASTA version will
        be produced in the same folder with the same name as the FASTQ file
        (but with .fasta instead of .fastq)
    '''
    def gene_calling(self, file, output, assembled = True, error_model = 'illumina_10'):
        bashCommand = 'run_FragGeneScan.pl -genome='
        if assembled:
            bashCommand += file + ' -out=' + output + '/fgs -complete=1 -train=./complete'
        else:
            mtools.fastq2fasta(file, file.replace('fastq', 'fasta'))            # fraggenescan only accepts FASTA input
            bashCommand += (file.replace('fastq', 'fasta') + ' -out=' + output +
                            '/fgs -complete=0 -train=./' + error_model)
        bashCommand += ' -thread=' + self.threads
        mtools.run_command(bashCommand)
    
    def annotation(self, max_target_seqs = '50'):
        diamond = DIAMOND(threads = self.threads,
                          db = self.db,
                          out = self.out_dir + '/aligned.blast',
                          query = self.out_dir + '/fgs.faa',
                          un = self.out_dir + '/unaligned.fasta',
                          unal = '1',
                          max_target_seqs = max_target_seqs)
        
        if self.db[-6:] == '.fasta':
            print('FASTA database was inputed')
            if not os.path.isfile(self.db.replace('fasta','dmnd')):
                print('DMND database not found. Generating a new one')
                diamond.set_database(self.db, self.db.replace('fasta','dmnd'))
            else:
                print('DMND database was found. Using it')
        elif self.db[-5:] != '.dmnd':
            print('Database must either be a FASTA (.fasta) or a DMND (.dmnd) file')
            exit()
        diamond.db = self.db.split('.dmnd')[0] if '.dmnd' in self.db else self.db.split('.fasta')[0]
        
        diamond.run()
        
    '''
    Input: 
        blast: name of the file from the DIAMOND annotation
        uniprotinfo: name of the uniprot information file
        cog_blast: name of the COG annotation file from self.annotate_cogs
        fun: name of the fun.txt file
        split_pathways: boolean, if MOSCA should split the information from the
        Pathway column provided by UniProt mapping
        readcounts_matrix: name of the file containing protein expression
    Output: 
        pandas.DataFrame with integration of information from UniProt and COGs
    '''
    def join_reports(self, blast, uniprotinfo, cog_blast, out_dir,
                     fun = 'MOSCA/Databases/COG/fun.txt', split_pathways = False):
        result = mtools.parse_blast(blast)
        result.index = result.qseqid
        result = result[result.sseqid != '*']
        result['Entry'] = [ide.split('|')[1] for ide in result.sseqid]
        uniprotinfo = pd.read_csv(uniprotinfo, sep= '\t', low_memory=False).drop_duplicates()
        if split_pathways:                                                      # TODO - the reorganization of pathways incurs in extra lines for same IDs. Find workaround
            print('Reorganizing pathways information.')
            funcdf = uniprotinfo[uniprotinfo.Pathway.notnull()][['Entry','Pathway']]
            funcdf.Pathway = funcdf.Pathway.apply(self.split)
            funcdf = self.using_repeat(funcdf)
            pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) 
                                        for path in funcdf.Pathway], index = funcdf.index)
            pathways.columns = ['Superpathway','Pathway','Subpathway']
            del funcdf['Pathway']; del uniprotinfo['Pathway']
            funcdf = pd.concat([pathways, funcdf], axis = 1)
            uniprotinfo = pd.merge(uniprotinfo, funcdf, on = ['Entry'], how = 'left')
        result = pd.merge(result, uniprotinfo, on = ['Entry'], how = 'left')
        cog_blast = pd.read_csv(cog_blast, sep = '\t')
        result = pd.merge(result, cog_blast, on = ['qseqid'], how = 'left')
        
        result.to_csv(out_dir + '/full_analysis_results.tsv', sep = '\t', index = False)
        
        print('Defining best consensus COG for each UniProt ID.')
        cogs_df = pd.DataFrame()
        tqdm.pandas()
        cogs_df['cog'] = result.groupby('Entry')['cog'].progress_apply(lambda x:
            x.value_counts().index[0] if len(x.value_counts().index) > 0 else np.nan)
        cogs_relation = result[['COG general functional category','COG functional category',
                            'COG protein description','cog']]
        cogs_df['Entry'] = cogs_df.index
        cogs_relation = cogs_relation.drop_duplicates()
        cogs_df = cogs_df[cogs_df['cog'].notnull()]
        cogs_df = pd.merge(cogs_df, cogs_relation, on = 'cog', how = 'inner')
        result.drop(['COG general functional category', 'COG functional category',
              'COG protein description', 'cog'], axis=1, inplace=True)
        result = pd.merge(result, cogs_df, on = 'Entry', how = 'outer')
        result = result[['Entry','Taxonomic lineage (SUPERKINGDOM)',
              'Taxonomic lineage (PHYLUM)','Taxonomic lineage (CLASS)',
              'Taxonomic lineage (ORDER)','Taxonomic lineage (FAMILY)',
              'Taxonomic lineage (GENUS)','Taxonomic lineage (SPECIES)',
              'EC number', 'Ensembl transcript','Function [CC]', 'Gene names',
              'Gene ontology (GO)', 'Keywords','Pathway', 'Protein existence',
              'Protein families', 'Protein names','Cross-reference (BRENDA)', 
              'Cross-reference (BioCyc)','Cross-reference (CDD)', 
              'Cross-reference (InterPro)','Cross-reference (KEGG)', 
              'Cross-reference (KO)','Cross-reference (Pfam)', 
              'Cross-reference (Reactome)','Cross-reference (RefSeq)', 
              'Cross-reference (UniPathway)','Cross-reference (eggNOG)',
              'COG general functional category', 'COG functional category',
              'COG protein description','cog']]
        result.columns = ['Protein ID'] + result.columns.tolist()[1:]
        result = result.drop_duplicates()
        
        result.to_csv(out_dir + '/protein_analysis_results.tsv', sep = '\t', index = False)
        return result
    
    def run(self):
        self.gene_calling(self.file, self.out_dir, self.assembled)
        #self.annotation()                                                      # annotation has to be refined to retrieved better than hypothetical proteins
        self.annotation(max_target_seqs = '1')
                    
    def global_information(self):
        # Join blast reports
        if not os.path.isfile(self.out_dir + '/Annotation/aligned.blast'):
            mtools.run_command('cat ' + ' '.join(glob.glob(self.out_dir + '/Annotation/*/aligned.blast')), 
                           file = self.out_dir + '/Annotation/aligned.blast')
        
        # Retrieval of information from UniProt IDs
        mtools.run_command('python UPIMAPI/upimapi.py -i {0}/Annotation/aligned.blast -o {0}/Annotation/uniprotinfo.tsv --full-id --blast{1}{2}'.format(
                self.out_dir, ' -anncols {}'.format(self.columns[0]) if self.columns != [''] else '',     # if columns are set, they will be inputed
                ' -anndbs {}'.format(self.databases[0]) if self.databases != [''] else ''))               # if databases are set, they will be inputed
        
        # Join COG reports
        if not os.path.isfile(self.out_dir + '/Annotation/protein2cog.tsv'):
            mtools.run_command('cat ' + ' '.join(glob.glob(
                    self.out_dir + '/Annotation/*/protein2cog.tsv')), 
                file = self.out_dir + '/Annotation/protein2cog.tsv')
        
        # Integration of all reports - BLAST, UNIPROTINFO, COG
        joined = self.join_reports(self.out_dir + '/Annotation/aligned.blast', 
                               self.out_dir + '/Annotation/uniprot_info.tsv', 
                               self.out_dir + '/Annotation/protein2cog.tsv',
                               self.out_dir)
        
        # MG quantification for each MG name of each Sample
        for sample in self.sample2name.keys():
            for mg_name in self.sample2name[sample]:
                mtools.perform_alignment('{}/Assembly/{}/contigs.fasta'.format(self.out_dir, sample),
                        ['{}/Preprocess/Trimmomatic/quality_trimmed_{}_{}_paired.fq'.format(self.out_dir, mg_name, fr)
                        for fr in ['forward', 'reverse']], '{}/Annotation/{}/{}'.format(self.out_dir, sample, mg_name),
                        threads = self.threads)
                mtools.normalize_readcounts_by_size('{}/Annotation/{}/{}.readcounts'.format(
                        self.out_dir, sample, mg_name), assembler = self.assembler)
                joined = mtools.define_abundance(joined, readcounts = '{}/Annotation/{}/{}.readcounts'.format(self.out_dir, sample, mg_name), 
                                             blast = '{}/Annotation/{}/aligned.blast'.format(self.out_dir, sample),
                                             name = mg_name)
        
        joined.to_csv(self.out_dir + '/joined_information.tsv', sep = '\t', index = False)
        print('joined was written to ' + self.out_dir + '/joined_information.tsv')  
        
        # TODO - this EXCEL writing is not working
        '''
        if os.path.isfile(self.out_dir + '/joined_information.xlsx'):
            os.remove(self.out_dir + '/joined_information.xlsx')
        
        # Write EXCEL - at the time of testing, maximmum lines for Excel file were 1048575 per sheet, so must split data between sheets # TODO - this takes forever, see what's wrong here
        writer = pd.ExcelWriter(self.out_dir + '/joined_information.xlsx', engine='xlsxwriter')
        i = 0
        j = 1
        while i + 1000000 < len(joined):
            joined.iloc[i:(i + 1000000)].to_excel(writer, sheet_name='Sheet ' + str(j), index = False)
            j += 1
        joined.iloc[i:len(joined)].to_excel(writer, sheet_name='Sheet ' + str(j), index = False)
        print('joined was written to ' + self.out_dir + '/joined_information.xlsx')
        '''
        self.joined2kronas(joined, self.out_dir + '/Annotation/krona', self.mg_names)
        
        return joined

    '''
    UniProt regularly updates its databases. As we are working with MG here, many
    identifications will sometimes pertain to ORFs or pseudogenes that have been 
    wrongly predicted to code for proteins. It may also happen that the original
    authors decided to withdraw their published sequences.
    Input:
        uniprotinfo: name of UniProt info file
    Output:
        Rows in UniProt info file that lack a species identification will be 
        updated to include that new information
    '''
    def info_from_no_species(self, uniprotinfo):
        uniprotinfo = pd.read_csv(uniprotinfo, sep = '\t')
        missing_uniprotinfo = uniprotinfo[uniprotinfo['Taxonomic lineage (SPECIES)'].isnull()]
        ids = list(set([ide for ide in missing_uniprotinfo['Entry']]))
        new_uniprotinfo = self.get_uniprot_information(ids)
        for entry in new_uniprotinfo['Entry']:
            missing_uniprotinfo[missing_uniprotinfo['Entry'] == entry] = new_uniprotinfo[new_uniprotinfo['Entry'] == entry]
        new_uniprotinfo.to_csv(uniprotinfo, sep = '\t', index = False)

    '''
    Input:
        tsv: filename of TSV file to be inputed. Must have the format 
        value\tcategorie1\tcategorie2\t..., with no header
        output: filename of HTML krona plot to output
    Output:
        A krona plot will be created at output if it has been specified, else
        at tsv.replace('.tsv','.html')
    '''
    def create_krona_plot(self, tsv, output = None):
        if output is None:
            output = tsv.replace('.tsv','.html')
        mtools.run_command('perl Krona/KronaTools/scripts/ImportText.pl {} -o {}'.format(tsv, output))
        
    '''
    Input:
        tsv: filename of MOSCA result from analysis
        output: basename for krona plots
        mg_columns: names of columns with abundances from which to build krona plots
    Output:
    '''
    def joined2kronas(self, joined, output, mg_columns, taxonomy_columns = [
            'Taxonomic lineage (' + level + ')' for level in ['SUPERKINGDOM', 
                    'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']],
                      functional_columns = ['COG general functional category',
                                            'COG functional category',
                                            'COG protein description']):
        print('Representing taxonomic and functional data in Krona plots at ' + output)
        for name in mg_columns:
            partial = joined.groupby(taxonomy_columns)[name].sum().reset_index()
            partial = partial[[name] + taxonomy_columns]
            partial.to_csv('{}_{}_tax.tsv'.format(output, name), sep = '\t', 
                           index = False, header = False)
            self.create_krona_plot('{}_{}_tax.tsv'.format(output, name))
            
            partial = joined.groupby(functional_columns)[name].sum().reset_index()
            partial = partial[[name] + functional_columns]
            partial.to_csv('{}_{}_fun.tsv'.format(output, name), sep = '\t', 
                           index = False, header = False)
            self.create_krona_plot('{}_{}_fun.tsv'.format(output, name))
            
    '''
    Input:
        data: str - result from MOSCA analysis
        seqs_file: str - file to store FASTA sequences
        blast_file: str - file to store extra annotation
    Output:
        returns data with new protein names
    '''
    def further_annotation(self, data, temp_folder = 'temp',
                           dmnd_db = 'MOSCA/Databases/annotation_databases/uniprot.dmnd',
                           threads = '12', out_format = '6', all_info = True):
        '''
        noid_entries = data[(data['Protein names']=='Uncharacterized protein') & 
                       (data['COG general functional category']=='POORLY CHARACTERIZED')]['Entry']

        self.recursive_uniprot_fasta(temp_folder + '/temp.fasta', 
                                                   entries = noid_entries)
        
        mtools.run_command(('diamond blastp -q {} --db {} -p {} -f {} -o {}').format(
                temp_folder + '/temp.fasta', dmnd_db, threads, out_format, 
                temp_folder + '/temp_blast.tsv'))
        '''
        blast = mtools.parse_blast(temp_folder + '/temp_blast.tsv')
        '''
        noid_entries = [ide.split('|')[1] for ide in blast['sseqid']]
        
        uniprotinfo = self.get_uniprot_information(noid_entries)
        uniprotinfo.to_csv(temp_folder + '/up.info',index=False)
        '''
        uniprotinfo = pd.read_csv(temp_folder + '/up.info')
        uniprotinfo.drop_duplicates(inplace = True)
        
        relation = pd.DataFrame()
        relation['qseqid'] = [ide.split('|')[1] for ide in blast['qseqid']]
        relation['sseqid'] = [ide.split('|')[1] for ide in blast['sseqid']]
        relation = pd.merge(relation, uniprotinfo, left_on = 'sseqid', 
                            right_on = 'Entry', how = 'inner')
        relation = relation[(relation['Protein names'] != 'Deleted.') &
                            (relation['Protein names'] != 'Conserved protein') &
                            (~relation['Protein names'].str.contains('uncharacterized', case = False))]
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)']
        data['Alternative entry'] = np.nan; data['Alternative protein name'] = np.nan
        pbar = ProgressBar()
        print('Searching for new annotations in the blast matches.')
        print('Proteins queried for new annotation:', len(set(relation['qseqid'])))
        print('Total possible identifications:', str(len(relation)))
        for query in pbar(relation['qseqid']):
            partial = relation[relation['qseqid'] == query]
            if all_info:
                data.loc[data['Entry'] == query, ['Alternative entry']] = '; '.join(partial['sseqid'])
                data.loc[data['Entry'] == query, ['Alternative protein name']] = '; '.join(partial['Protein names'])
            else:
                if len(partial) > 0:
                    original_taxonomy = data.loc[data['Entry'] == query].iloc[0][tax_columns]
                    max_proximity = 0; entry = np.nan; protein_name = np.nan
                    for i in range(len(partial)):
                        tax_score = self.score_taxonomic_proximity(original_taxonomy,
                            partial.iloc[i][tax_columns])
                        if tax_score > max_proximity:
                            max_proximity = tax_score
                            entry = partial.iloc[i]['sseqid']
                            protein_name = partial.iloc[i]['Protein names']
                    data.loc[data['Entry'] == query, ['Alternative entry']] = entry
                    data.loc[data['Entry'] == query, ['Alternative protein name']] = protein_name
        return data

    def score_taxonomic_proximity(self, tax1, tax2):
        proximity = 0
        for level in range(len(tax1)):
            if tax1[level] != tax2[level]:
                return proximity
            else:
                proximity += 1
        return proximity