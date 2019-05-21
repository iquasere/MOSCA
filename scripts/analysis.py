# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package

By João Sequeira

Apr 2019
"""

'''
COG articles
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563889/
-They compared COG composition obtained from their organisms with COG composition of random organisms
    -it revealed a shift towards energy production (more COGs from those categories)

Significant differences
-first, all of the genes significant different through considering the duplicates
-after that, the genes of specific COGs

Tools for mappping KEGG(or not) metabolic pathways
https://cytoscape.org/what_is_cytoscape.html
    http://iapg2p.sourceforge.net/kgml-ed/Introduction.html
http://www.ra.cs.uni-tuebingen.de/software/KEGGtranslator/
https://pathways.embl.de/
http://www.g-language.org/data/marray/software/map2swf.cgi
https://biit.cs.ut.ee/kegganim/
'''

import pandas as pd
import numpy as np
from scipy import stats

class Analysis:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        self.joined = self.parse_joined_output(self.file)
    
    '''
    input:
        file: joined CSV output file from MOSCA's workflow
    output
        pd.DataFrame object
    '''
    def parse_joined_output(self, file):
        return pd.read_csv(file, sep = '\t')                                    # TODO - check if needed to name sheet to import first sheet
    
    def t_test(self, joined):
        pathways = ['One-carbon metabolism; methanogenesis from CO(2)',
                    'One-carbon metabolism; methanogenesis from acetate',
                    'Cofactor biosynthesis', 'lipid metabolism', 
                    'Metabolic intermediate biosynthesis']
        protein_names = ['ATP synthase']
        
        rows_of_interest = joined[]
        

mgmp[['Protein ID']+mgmpmg_columns].to_csv('MGMP/de_protein_name_mg.ssv',sep=' ')
mgmp[['Protein ID']+[col + '_spectra' for col in mgmpmg_columns]].to_csv('MGMP/de_protein_name_mp.ssv',sep=' ')

mgmt[['Protein ID']+mgmtmg_columns].to_csv('MOSCAfinal/de_protein_name_mg.ssv',sep=' ')
mgmt[['Protein ID']+mgmtmt_columns].to_csv('MOSCAfinal/de_protein_name_mt.ssv',sep=' ')


mgmp_cogs = mgmp[cog_columns + mgmpmg_columns + mgmpmp_columns]
mgmt_cogs = mgmt[cog_columns + mgmtmg_columns + mgmtmt_columns]

mgmp_cogs = mgmp_cogs.groupby(cog_columns)[mgmpmg_columns + spectra_columns].sum().reset_index()
mgmt_cogs = mgmt_cogs.groupby(cog_columns)[mgmtmg_columns + mgmtmt_columns].sum().reset_index()


mgmt_cogs.groupby(cog_columns[0])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl0_mg.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[1])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl1_mg.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[2])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl2_mg.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[3])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl3_mg.ssv',sep=' ',index=False)

mgmt_cogs.groupby(cog_columns[0])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl0_mt.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[1])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl1_mt.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[2])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl2_mt.ssv',sep=' ',index=False)
mgmt_cogs.groupby(cog_columns[3])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/de_cogs_lvl3_mt.ssv',sep=' ',index=False)


mgmp_cogs.groupby(cog_columns[0])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl0_mg.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[1])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl1_mg.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[2])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl2_mg.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[3])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl3_mg.ssv',sep=' ',index=False)

mgmp_cogs.groupby(cog_columns[0])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl0_mp.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[1])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl1_mp.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[2])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl2_mp.ssv',sep=' ',index=False)
mgmp_cogs.groupby(cog_columns[3])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/de_cogs_lvl3_mp.ssv',sep=' ',index=False)

mgmp_cogs.groupby(cog_columns)[mgmpmg_columns+mgmpmp_columns].sum().reset_index().to_excel('MGMP/krona_cog.xlsx',index=False)
mgmp_tax.groupby(tax_columns)[mgmpmg_columns+mgmpmp_columns].sum().reset_index().to_excel('MGMP/krona_tax.xlsx',index=False)

mgmt_cogs.groupby(cog_columns)[mgmtmg_columns+mgmtmt_columns].sum().reset_index().to_excel('MOSCAfinal/krona_cog.xlsx',index=False)
mgmt_tax.groupby(tax_columns)[mgmtmg_columns+mgmtmt_columns].sum().reset_index().to_excel('MOSCAfinal/krona_tax.xlsx',index=False)


mgmp_tax = mgmp[tax_columns + mgmpmg_columns + mgmpmp_columns]
mgmt_tax = mgmt[tax_columns + mgmtmg_columns + mgmtmt_columns]

mgmt_tax.groupby(tax_columns[0])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl0_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[1])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl1_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[2])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl2_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[3])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl3_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[4])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl4_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[5])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl5_mg.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[6])[mgmtmg_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl6_mg.ssv',sep=' ',index=False)

mgmt_tax.groupby(tax_columns[0])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl0_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[1])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl1_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[2])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl2_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[3])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl3_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[4])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl4_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[5])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl5_mt.ssv',sep=' ',index=False)
mgmt_tax.groupby(tax_columns[6])[mgmtmt_columns].sum().reset_index().to_csv('MOSCAfinal/tax_lvl6_mt.ssv',sep=' ',index=False)


mgmp_tax.groupby(tax_columns[0])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl0_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[1])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl1_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[2])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl2_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[3])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl3_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[4])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl4_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[5])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl5_mg.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[6])[mgmpmg_columns].sum().reset_index().to_csv('MGMP/tax_lvl6_mg.ssv',sep=' ',index=False)

mgmp_tax.groupby(tax_columns[0])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl0_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[1])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl1_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[2])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl2_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[3])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl3_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[4])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl4_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[5])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl5_mp.ssv',sep=' ',index=False)
mgmp_tax.groupby(tax_columns[6])[mgmpmp_columns].sum().reset_index().to_csv('MGMP/tax_lvl6_mp.ssv',sep=' ',index=False)



