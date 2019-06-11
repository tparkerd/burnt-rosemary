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
  """Validates input file for phenotype data

  This function validates that the contents of a file to contain phenotype data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
  df = pd.read_csv(filepath)
  nrows, ncols = df.shape
  nrows += 1 # include the header in the row count

  schema_columns = [
    Column('Pedigree', [
      IsDistinctValidation()
    ])
  ]

  for n in range(1, ncols):
    schema_columns.append(
      Column(df.columns[n], [
        # NOTE(tparker): This may not always be true. If there any phenotypes that
        # are listed as categories or strings, then this would fail
        # Find out all the possible phenotype values. It may be difficult to
        # validate input data without a user-provided dtype list
        CanConvertValidation(float)
      ])
    )

  schema = Schema(schema_columns)
  err = schema.validate(df)

  if err:
    for e in err:
      logging.error(e)
      raise Exception(e)

def validate_line(conn, args, filepath):
  """Validates input file for line data

  This function validates that the contents of a file to contain line data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
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
      raise Exception(e)

def validate_genotype(conn, args, filepath):
  """Validates input file for genotype data

  This function validates that the contents of a file to contain genotype data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
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
      raise Exception(e)

def validate_variant(conn, args, filepath):
  """Validates input file for variant data

  This function validates that the contents of a file to contain variant data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
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
      raise Exception(e)

def validate_kinship(conn, args, filepath):
  """Validates input file for kinship data

  This function validates that the contents of a file to contain kinship data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
  df = pd.read_csv(filepath)
  nrows, ncols = df.shape
  df.rename(columns = {"Unnamed: 0": "line_name"}, inplace=True) # since column name is blank by default, rename it for later reference
  nrows += 1 # include the header row in the count
  logging.debug(f"Dimensions of kinship matrix: <{nrows}, {ncols}>")

  schema_columns = [
    Column('line_name', [
      IsDistinctValidation()
    ])
  ]

  for n in range(1, ncols):
    schema_columns.append(Column(df.columns[n], [
      CanConvertValidation(float)
    ]))

  schema = Schema(schema_columns)

  err = schema.validate(df)

  if err:
    for e in err:
      logging.error(e)
      raise Exception(e)
  
def validate_population_structure(conn, args, filepath):
  """Validates input file for population structure data

  This function validates that the contents of a file to contain population
  structure data. If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
  df = pd.read_csv(filepath)
  nrows, ncols = df.shape
  nrows += 1 # include the header rows in the count
  logging.debug(f'Population structure columns: {df.columns}')
  logging.debug(f"Population structure dimensions: <{nrows}, {ncols}>")


  schema_columns = [
    Column('Pedigree', [
      IsDistinctValidation()
    ])
  ]

  for n in range(1, ncols):
    schema_columns.append(Column(df.columns[n], [
      CanConvertValidation(float)
    ]))

  schema = Schema(schema_columns)
  err = schema.validate(df)

  if err:
    for e in err:
      logging.error(e)
      raise Exception(e)

def validate_runs(conn, args, filepath):
  """Validates input file for GWAS run data

  This function validates that the contents of a file to contain GWAS run data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
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
      raise Exception(e)

def validate_results(conn, args, filepath):
  """Validates input file for GWAS result data

  This function validates that the contents of a file to contain GWAS result data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    filepath (str): location of input file
  
  """
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
      raise Exception(e)


