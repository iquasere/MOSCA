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