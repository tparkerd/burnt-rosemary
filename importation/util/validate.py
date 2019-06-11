"""Validation methods for importation"""
import csv
import errno
import itertools
import logging
import os
import re

import numpy as np
import pandas as pd
import psycopg2
from pandas_schema import Column, Schema
from pandas_schema.validation import (CanConvertValidation,
                                      CustomSeriesValidation,
                                      InRangeValidation, IsDistinctValidation,
                                      MatchesPatternValidation)

import importation.util.insert as insert
from importation.util.dbconnect import config, connect
from importation.util.models import (chromosome, genotype, genotype_version,
                                     growout, growout_type, gwas_algorithm,
                                     gwas_result, gwas_run, imputation_method,
                                     kinship, kinship_algorithm, line,
                                     location, phenotype, population,
                                     population_structure,
                                     population_structure_algorithm, species,
                                     trait, variant)


def file_exists(conn, args, filepath):
  if not os.path.isfile(filepath):
    return True
  else:
    logging.error("%s: %s", os.strerror(errno.ENOENT), filepath)
    return False

# Check that a growout location exists in the database
def location_exists(conn, args, location):
  pass

# PHENOTYPE
# Check that the phenotype file has required fields
# Check if it is a combined growout phenotype file
def validate_phenotype(conn, args, filepath):
  pass

# GENOTYPE
# Check that the genotype file has required fields or at least has the correct
# file format (012)
def validate_line(conn, args, filepath):
  schema = Schema([
    Column('line_name', [
      IsDistinctValidation()
    ])
  ])

  df = pd.read_csv(filepath, header=None)

  if len(df.columns) != 1:
    raise Exception(f"Invalid file format. Excepted 1 column found {len(df.columns)} columns. This file should be a single column of each line. Each entry should be distinct.")
  
  df.columns = [ 'line_name' ]
  err = schema.validate(df)

  if err:
    for e in err:
      logging.error(e)

def validate_genotype(conn, args, filepath):
  # Allow for users to skip this validation step because it is time consuming
  if args.skip_genotype_validation is True:
    return

  
  schema_columns = [
    Column('row_number', [
      CanConvertValidation(int) &
      IsDistinctValidation()
    ])
  ]

  # Get the number of lines from the .pos counterpart file
  pos_filepath = '.'.join([filepath, 'pos'])
  if not os.path.exists(pos_filepath):
    raise FileNotFoundError(f"Count not locate the position counterpart file for {filepath}")
  nPositions = len(pd.read_csv(pos_filepath, header=None).index)

  for n in range(nPositions):
    schema_columns.append(
      Column(f'pos_{n}', [
        CanConvertValidation(int) &
        CustomSeriesValidation(lambda x: x.int in [-1,0,1,2], 'Incorrectly coded value.')
      ])
    )

  schema = Schema(schema_columns)

  df = pd.read_csv(filepath, sep='\t', header=None)

  err = schema.validate(df)
  
  if err:
    for e in err:
      logging.error(e)

def validate_variant(conn, args, filepath):
  schema = Schema([
    Column('chr', [
      CanConvertValidation(int)
    ]),
    Column('pos', [
      CanConvertValidation(int)
    ])
  ])

  df = pd.read_csv(filepath, sep='\t', header=None)

  if len(df.columns) != 2:
    raise Exception(f"Invalid file format. Excepted 2 columns, found {len(df.columns)} columns. Columns should consist of chromsome number and SNP position. Filepath: {filepath}")

  df.columns = ['chr', 'pos']
  err = schema.validate(df)

  if err:
    for e in err:
      logging.error(e)

def validate_kinship(conn, args, filepath):
  df = pd.read_csv(filepath)
  nrows, ncols = df.size
  logging.info("%s, %s", nrows, ncols)
  
def validate_population_structure(conn, args, filepath):
  pass

# GWAS RUNS
# Check that the GWAS runs file contains the required fields
# trait
def validate_runs(conn, args, filepath):
  pass

# GWAS RESULTS
# Check that the GWAS results file contains the required fields
# SNP, nSNPs, chr, pos
def validate_results(conn, args, filepath):

  df = pd.read_csv(filepath)
  # For each column, add it to the schema, and then for known ones, add the 
  # schema validation. Use fuzzy comparisons when possible
  schema_columns = []
  for col in df.columns:
    validators = []
    if re.match("(SNP)|(chr)|(pos)|(nSNPs)", col, re.IGNORECASE):
      validators.append(CanConvertValidation(int))
    # Look for any of the p-values and make sure that they can be cast as a float
    if re.match("((null)?pval(ue)?)", col, re.IGNORECASE):
      validators.append(CanConvertValidation(float))
    
    schema_columns.append(Column(col, validators))
  schema = Schema(schema_columns)

  err = schema.validate(df)
  if err:
    for e in err:
      logging.error(e)


