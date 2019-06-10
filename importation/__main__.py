#!/usr/bin/env python
"""Importation script that handles retrieving user-input for initializing data"""

import argparse
import os
import sys
import shutil
import datetime
import json
import importation.util.importutility as iu
import importation.util.find as find
from importation.util.dbconnect import connect

from pprint import pprint

def process(args, conn):
  """Method to preprocess dataset to parse out the necessary """
  # try:
  #   with open(args.meta, 'r') as mfp:
  #     m = json.load(mfp)

  #     pprint(m)
  # except:
  #   raise

  # NOTE: Example of connecting to the database to pull a bit of info
  # This can also be returned from the database as part of the sql statement
  # like the ones that Molly wrote
  # # Okay, SSH tunneled and it works just fine
  # # I mean, I could implement my own SSH tunneling like DBeaver
  # cur = conn.cursor()
  # sql = """SELECT * FROM chromosome;"""
  # cur.execute(sql)
  # rows = cur.fetchall()
  # for row in rows:
  #   pprint(row)
  # cur.close()

  pass


def init(args):
  """Initialization of a dataset. This function will initialize a configuration
  file for the purpose of preprocessing the data and offering a pre-populated
  list of the files to be used."""

  print("""This utility will walk you through creating a config.json file.

Enter any known values. Any changes can be made post-initialization by
altering the generated JSON file.
  
Press ^C at any time to quit.""")

  # Initialize configuration and store in dictory
  # Directory will be converted to JSON object and saved as config.json
  config = {}
  # Species & Population
  config["species"] = {}
  config["species"]["shortname"] = input(f"species shortname: ") or ""
  config["species"]["binomial name"] = input(f"species binomial name: ") or ""
  config["population"] = input(f"population: ") or ""
  # Growout(s)
  config["growout"] = []
  growout_index = 0
  while (True):
    config["growout"].append({})
    config["growout"][growout_index]["name"] = input(f"growout name #{growout_index + 1}: ") or ""
    suggested_year = iu.growout_name_to_year(config["growout"][growout_index]["name"])
    if suggested_year is not None:
      config["growout"][growout_index]["year"] = input(f"growout year #{growout_index + 1} ({suggested_year}): ") or str(suggested_year)
    else:
      config["growout"][growout_index]["year"] = input(f"growout year #{growout_index + 1}: ") or ""
      
    config["growout"][growout_index]["type"] = input(f"growout type #{growout_index + 1}: ") or ""
    config["growout"][growout_index]["location"] = {}
    config["growout"][growout_index]["location"]["code"] = input(
        f"growout location code #{growout_index + 1}: ") or ""
    growout_index += 1
    stdin = input(f"Would you like to add another growout? [y/N]: ").lower() or 'n'
    if stdin == 'n' or stdin == "no":
      break

  # Check if there is already a configuration file
  # If there is one, prompt user to overwrite it
  try:
    if os.path.exists(args.config):
      stdin = input(
          f"File '{args.config}' already exists. Do you want to overwrite it? [y/N]: ").lower() or "n"
      if stdin == 'n' or stdin == "no":
        print("Initialization aborted.")
        sys.exit(0)
    print(f"About to write to {args.config}:\n")
    pprint(config)
    stdin = input(f"\nIs this okay? [Y/n]: ").lower() or "y"
    if stdin == 'n' or stdin == "no":
      print("Initialization aborted.")
      sys.exit(0)

    if not args.debug:
      with open(args.config, 'w') as ofp:
        json.dump(config, ofp, indent=2)
  except:
    raise


# ================== Design Stage ==================
def design(args):
  """Validates JSON configuration file and checks if experiement can be created
  by checking against database."""

  cwd = os.getcwd()
  try:
    with open(args.config, 'r', encoding='utf-8') as configfp:
      config = json.load(configfp)
  except:
    raise
  
  conn = connect()

  if args.debug:
    pprint(conn)

  # Species & Poputation
  # Check if the species is already in the database
  species = find.find_species(conn, config["species"]["shortname"])
  if args.debug:
    print(f"Species: {config['species']['shortname']} -> {species}")

  # Growouts
  for go in config["growout"]:
    go["id"]= find.find_growout(conn, go["name"])

    location = {}
    if "code" in go["location"]:
      location_id = find.find_location(conn, go["location"]["code"])
      if location_id is not None:
        # Location exists, no further action required
        # Go to the next growout
        continue
      else:
        location = iu.create_location(go["location"]["code"])
    # Unknown location: location not in database OR location code not provided
    # If the code is provided but not in database, create a new entry for the
    # database
    if location_id is None:
      location = iu.create_location(go["location"]["code"])

    # Else, code is not provided, so we need to continually request that it be
    # defined and then entered into the database


    while True:
      code = input(f"Define unique location code: ").strip()
      if code is None or code == "":
        print("Location code cannot be blank. Please enter an alphanumeric string.")
        continue
      elif find.find_location(conn, code) is not None:
        print(f"Location '{code}' is already in use. Please enter a different code.")
        continue
      else:
        break
       
      go["location"]["id"] = find.find_location(conn, go["location"]["code"])
        

    
    pprint(go)
    
  # Location

  # Year
  # Type


# ================== Collection Stage ==================
def collect(args):
  """Collects additional information on experiment and how the GWAS will be executed
  
  This handles importing the phenotypic data
  
  """
  pass


# ================== Results Stage ==================
def result(args):
  """Collects the metadata and resultant dataset of a GWAS run

  This handles collecting user input about a GWAS run and then associates it with
  a result set and then imports it.
  """

  pass




def parseOptions():
  """
  Function to parse user-provided options from terminal
  """
  description='Importation script that handles retrieving user-input for initializing data.'
  usage="python import.py -i ./data_file_directory [-o ./output_directory]"
  default_config_location = os.path.join(os.getcwd(), 'config.json')
  parser = argparse.ArgumentParser(description=description, usage=usage)
  parser.add_argument("--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("-o", "--outdir", default = f"output_{datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')}", help="Path of output directory")
  parser.add_argument("--debug", action="store_true", help="Enables --verbose and disables writes to disk")
  parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0-alpha')
  parser.add_argument("-i", "--input_directory", default="config.json", help="Input directory of data files. See documentation on required file structure hierarchy.")
  parser.add_argument("-c", "--config", default=default_config_location, help="JSON format file that contains metadata about the dataset to import.")
  parser.add_argument("--init", action="store_true", help="Scans directory for all files required to perform importation. Generates suggested JSON file used for importation")
  parser.add_argument("--design", action="store_true", help="Takes JSON file and validates if experiment can be imported")
  parser.add_argument("--collect", action="store_true", help="Gathers all the additional information in preperation for running GWAS on the data")
  parser.add_argument("--result", action="store_true", help="Collects last bit of user input to import the GWAS run and results information")
  args = parser.parse_args()
  if args.debug is True:
    args.verbose = True
    args.write = False

  if args.debug:
    pprint(args)
  
  return args

if __name__ == "__main__":
  args = parseOptions()
  try:
    if args.init:
      init(args)
    if args.design:
      design(args)
    if args.collect:
      collect(args)
    if args.result:
      result(args)
    # conn = connect()
    # process(args, conn)
  except:
    raise


# Dependency rundown
# Independent values -- require extraction or user-input
#   NAME                           ABBR     DEP(S)
#   Location                       L
#   Species                        S
#   Traits                         T
#   Growout_type                   GT
#   GWAS_algorithm                 GA
#   Imputation_method              IM
#   Kinship_algorithm              KA
#   Population_structure_algortihm PSA
#
# One Independent Dependency
#   Population                     P        S
#   Chromosome                     C        S
#   Kinship                        K        KA
#   Population_structure           PS       PSA
#
# One Dependent Dependency
#   Line                           LI        P
#
# Two Dependencies (Mixed only)
#   Genotype_version               GV       L, P
#   Variant                        V        S, C
#   Phenotype                      PH       L, T
#
# Three Dependencies (Mixed only)
#   Genotype                       G        L, C, GV
#   Growout                        GO       LOC, P, GT
#
# Many Dependencies (Mixed only)
#   GWAS_run                       GRUN     T, GA, GV, IM, K, PS
#   GWAS_result                    GRES     C, GRUN

# ============= Independent Values =============
# Get all of the independent values
# NOTE: this includes finding them if the value already exists in the database
# Location
#   req. user input: country
#   These values can be extracted from phenotype files
#   Their filenames include the location code
# NOTE(timp): 'Combined' location exists, but unclear what it represents
#

# Species
#   req. user input: shortname, binomial
# Traits
#   req. user input: none
# Growout_type
#   req. user input: name
# GWAS_Algorithm
#   req. user input: name
# Imputation_method
#   req. user input: name
# Kinship_algorithm
#   req. user input: name
# Population_structure_algorithm
#   req. user input: name

# ============= Single 'Independency' =============
# Population
#   req. user input: name
# Chromosome
#   req. user input: ??? (uncertain if there is a naming convention)
# Kinship
#   req. user input: none (filepath generated by this script/target host)
# Population_structure
#   req. user input: none (filepath generated by this script/target host)

# ============= Single Dependency =============
# Line
#   req. user input: none (lines can be pulled from 012 files)

# ============= Two Dependencies =============
# Genotype_version
#   req. user input: ??? (does this require a lookup of the line?)
# Variant
#   req. user input: none (variants can be pulled from 012.pos file)
# Phenotype
#   req. user input: none (phenotype measurement value can be pulled from parsed pheno files)
#   The parsed pheno files originated from CSV transform on 5.mergedWeightNorm...longFormat.csv

# ============= Three Dependencies =============
# Genotype
#   req. user input: none (genotypes can be pulled from .012 and .012.indv files)
#   This will happen for each chromosome (and possibly scaffold) specified prior
# Growout
# NOTE(timp): Consider adding a new constraint to the growout to guarantee the uniqueness
#             of (population, location, year, and growout_type)
#             Should this mean that the name should *not* be unique?
#

# ============= Many Dependencies =============
# GWAS_run
#   req. user input: missing SNP cutoff value, missing line cutoff value, minor allele cutoff value
# GWAS_result
# NOTE(timp): For some reason the minor allele cutoff value was placed in a different spot than that in the GWAS_run table
#   req. user input: missing SNP cutoff value, missing line cutoff value, minor allele cutoff value

# Sequence of Events
# 1. Define an experiment/growout
#   - Location
#     - country name, code
#   - Species
#     - name
#   - Population
#   - Growout type
#   -
