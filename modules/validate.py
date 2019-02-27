"""Validation methods for importation"""
import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
from dbconnect import config, connect
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result
