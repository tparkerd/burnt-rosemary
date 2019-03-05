#!/bin/python
"""Hard-coded importation script for Setaria Diverity Panel dataset.
The goal of this is to just get the data into the database and have a better grasp on how to generalize this
"""
import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
import find
from dbconnect import config, connect
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result


# ========== Experiment Design ==========
# What the database needs in order to create an 'experiment' is the follow
# Species: maize (Zea mays)
# Population: Maize282
# Chromosome: 10 (just the number and a way to generate its unique name)
# Line: 282set_B73 (B73) -- taken from file if possible
# Genotype Version: B73 RefGen_v4_AGPv4_Maize282 (the reference genome)
# Growout, Type, and Location: 
#       Location: code, city, state, country
#                 "PU", "West Lafayette", "Indiana", "United States"
#       Type: "field", "phenotyper", etc.
#       Growout: name, population ID, location ID, year, type
#                "PU09", maize282popID, PUlocID, 2009, fieldGrowoutTypeID
# Traits (planned phenotypes/traits to measure)


# ========== Pipeline Design ==========
# GWAS Algorithm: "MLMM", "EMMAx", "GAPIT", "FarmCPU"
# Imputation Method: "impute to major allele"
# Kinship Algorithm: "loiselle"
# Population Structure Algorithm: "Eigenstrat"


# ========== Experiment Collection ==========
# Phenotype (external source?)
#       This needs to be standardized to a .pheno filetype.
#       For now, it is the longFormat for the Maize282 datset
        # 5.mergedWeightNorm.LM.rankAvg.longFormat.csv, but for Setaria will be 
# Genotype (VCF output)
#       
# Variant (VCF output)

# ========== Pipeline Collection ==========
# Kinship
# Population Structure

# ========== Results ==========
# GWAS Run
# GWAS Results
