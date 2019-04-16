# -*- coding: utf-8 -*-
"""
MOSCA's Analysis package

By João Sequeira

Apr 2019
"""

'''
COG articles
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5563889/
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
        file: joined output file from MOSCA's workflow
    output
        pd.DataFrame object
    '''
    def parse_joined_output(self, file):
        return pd.read_excel(file)                                            # TODO - check if needed to name sheet to import first sheet
    
    def t_test(self, joined):
        pathways = ['One-carbon metabolism; methanogenesis from CO(2)',
                    'One-carbon metabolism; methanogenesis from acetate',
                    'Cofactor biosynthesis', 'lipid metabolism', 
                    'Metabolic intermediate biosynthesis']
        protein_names = ['ATP synthase']
        
        rows_of_interest = joined[]