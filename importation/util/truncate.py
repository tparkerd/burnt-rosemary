"""Truncate the QA database"""

import logging
import time

import numpy as np
import pandas as pd
import psycopg2 as pg
from tqdm import tqdm

from importation.util.dbconnect import connect


def truncate(args):
  """This function resets by truncating all tables in the QA server
  
  Tables affected:
    growout_type
    line
    location
    growout
    variant
    genotype_version
    genotype
    kinship_algorithm
    phenotype
    trait
    population_structure_algorithm
    gwas_algorithm
    kinship
    population_structure
    gwas_run
    gwas_result
    tissue
    chromosome
    species
    population
    imputation_method
  
  """

  def empty(conn, args, table_name):
    """Remove all entries in a table with truncation

    This function removes all data in a table

    Args:
      conn (psycopg2.extensions.connection): psycopg2 connection
      args (ArgumentParser namespace): user-defined arguments
      table_name (string): name of table in database
    """
    cur = conn.cursor()

    # See if data has already been inserted, and if so, return it
    SQL = """TRUNCATE TABLE %s"""
    args_tuple = table_name
    try:
      cur.execute(SQL, args_tuple)    
    except:
      raise

  conn = connect(args)
  
  tables = [ 'growout_type',
             'line',
             'location',
             'growout',
             'variant',
             'genotype_version',
             'genotype',
             'kinship_algorithm',
             'phenotype',
             'trait',
             'population_structure_algorithm',
             'gwas_algorithm',
             'kinship',
             'population_structure',
             'gwas_run',
             'gwas_result',
             'tissue',
             'chromosome',
             'species',
             'population',
             'imputation_method' ]

  for table in tables:
    empty(conn, args, table)