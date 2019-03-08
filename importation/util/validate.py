"""Validation methods for importation"""
import pandas as pd
import numpy as np
import psycopg2
import csv
import util.insert as insert
from util.dbconnect import config, connect
from util.models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result

def location_exists(conn, location):
    pass