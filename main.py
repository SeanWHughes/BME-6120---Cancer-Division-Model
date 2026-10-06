# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 17:18:28 2026

@author: Sean
"""

#%%     IMPORT MODULES

# Standard modules
from dataclasses import dataclass, field

#%%     CLASSES

class Cell:
    def __init__(self):
        """
        Constructor method used to initialize new cell object instances.
        """
            self.celltype = 
            self.mut_rate = 0.05
            self.grow_rate = 
            self.death_rate = 0.05
            self.plasticity = 0
            self.genes = {}
            
            # Gene assignments
            genes['BRCA1'] = Gene('BRCA1', 'driver')
            genes['KRAS'] = Gene('KRAS', 'driver')  
            genes['ALB'] = Gene('ALB', 'neutral')
            genes['PAH'] = Gene('PAH', 'neutral')
            genes['APOA1'] = Gene('APOA1', 'neutral')
            
            genes['TTR'] = Gene('TTR', 'neutral')
            genes['TF'] = Gene('TF', 'neutral')
            genes['SERPINA1'] = Gene('SERPINA1', 'neutral')
            genes['CFTR'] = Gene('CFTR', 'neutral')
            genes['SLC2A2'] = Gene('SLC2A2', 'neutral')
            
            genes['SLC5A1'] = Gene('SLC5A1', 'neutral')
            genes['KCNJ10'] = Gene('KCNJ10', 'neutral')
            genes['AQP1'] = Gene('AQP1', 'neutral')
            genes['CA2'] = Gene('CA2', 'neutral')
            genes['ALDH2'] = Gene('ALDH2', 'neutral')
            
            genes['ADH1B'] = Gene('ADH1B', 'neutral')
            genes['UGT1A1'] = Gene('UGT1A1', 'neutral')
            genes['CYP2C19'] = Gene('CYP2C19', 'neutral')
            genes['CYP2D6'] = Gene('CYP2D6', 'neutral')
            genes['OTC'] = Gene('OTC', 'neutral')

class Gene:
    def __init__(self, gene_name: str, gene_purpose: str):
        """
        Constructor method used to initialize new gene object instances.
        """
        
        self.gene_name = gene_name
        self.gene_purpose = gene_purpose
        self.mutated = False

class Tissue:
    def __init__(self, tissue_name: str):
        """
        Constructor method used to initialize new gene object instances.
        """
        
        self.tissue_name = tissue_name
        
        
