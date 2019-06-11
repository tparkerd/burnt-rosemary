#!/usr/bin/env python
"""Importation script that handles retrieving user-input for initializing data"""

import argparse
import datetime
import errno
import json
import logging
import os
import shutil
import sys
from pprint import pprint
from string import Template

import importation.util.find as find
import importation.util.importutility as iu
from importation.util.dbconnect import connect
from importation.util.validate import (validate_genotype, validate_kinship,
                                       validate_line, validate_phenotype,
                                       validate_population_structure,
                                       validate_results, validate_runs,
                                       validate_variant)


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
  
  conn = connect(args)

  if args.debug:
    pprint(conn)

  # Species & Poputation
  # Check if the species is already in the database
  species = find.find_species(conn, config["species"]["shortname"])
  if args.debug:
    print(f"Species: {config['species']['shortname']} -> {species}")

  # Growouts
  for go in config["growout"]:
    go["id"]= find.find_growout(conn, args, go["name"])

    location = {}
    if "code" in go["location"]:
      location_id = find.find_location(conn, args, go["location"]["code"])
      if location_id is not None:
        # Location exists, no further action required
        # Go to the next growout
        continue
      else:
        location = iu.create_location(go["location"]["code"], args)
    # Unknown location: location not in database OR location code not provided
    # If the code is provided but not in database, create a new entry for the
    # database
    if location_id is None:
      location = iu.create_location(go["location"]["code"], args)

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


def validate(args):
  """Validate specified input files for import into GWAS database

  This function validates that the contents of a file to contain correct data.
  If an error is encountered, throw an exception.

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
  
  """
  logging.info("Collecting and identifying input files.")
  try:
    conn = connect(args)
  except:
    raise

  # Input file preprocessing and validation
  try:
    if not os.path.isfile(args.filename):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.filename)
    else:
      with open(args.filename) as f:
        dp = json.load(f) # data parameters

    # Verify that all necessary values were provided, assuming a complete dataset
      expected_fields = [ "species_shortname",
                        "species_binomial_name",
                        "species_subspecies",
                        "species_variety",
                        "population_name",
                        "number_of_chromosomes",
                        "genotype_version_assembly_name",
                        "genotype_version_annotation_name",
                        "reference_genome_line_name",
                        "gwas_algorithm_name",
                        "imputation_method_name",
                        "kinship_algortihm_name",
                        "population_structure_algorithm_name",
                        "kinship_filename",
                        "population_structure_filename",
                        "gwas_run_filename",
                        "gwas_results_filename",
                        "missing_SNP_cutoff_value",
                        "missing_line_cutoff_value",
                        "minor_allele_frequency_cutoff_value",
                        "phenotype_filename"
                      ]

    missing_keys = []
    for k in expected_fields:
      if k not in dp:
        missing_keys.append(k)
    if missing_keys:
      raise KeyError(f'The following keys are required. Please include them in your JSON configuration: {missing_keys}')

    # Check for all required fields
    required_fields = [ "species_shortname",
                        "species_binomial_name",
                        "population_name",
                        "number_of_chromosomes",
                        "genotype_version_assembly_name",
                        "genotype_version_annotation_name",
                        "reference_genome_line_name",
                        "gwas_algorithm_name",
                        "imputation_method_name",
                        "kinship_algortihm_name",
                        "population_structure_algorithm_name",
                        "kinship_filename",
                        "population_structure_filename",
                        "gwas_run_filename",
                        "gwas_results_filename",
                        "missing_SNP_cutoff_value",
                        "missing_line_cutoff_value",
                        "minor_allele_frequency_cutoff_value",
                        "phenotype_filename"
                      ]

    empty_fields = []
    for rf in required_fields:
      if not dp[rf]:
        empty_fields.append(rf)
    if empty_fields:
      raise KeyError(f'The following keys must be defined. Empty strings are not permitted. Please modify your JSON configuration: {empty_fields}')

    logging.info('Configuration file is valid. Verifying that all files exist.')

    # Track all the files to check for existance
    locations = []
    filepath_template = Template('${cwd}/${filename}')

    # Verify that all files exist
    # Lines
    lines_filename = Template('${chr}_${shortname}.012.indv')
    # Genotype
    genotype_filename = Template('${chr}_${shortname}.012')
    # Variants
    variants_filename = Template('${chr}_${shortname}.012.pos')

    for c in range(1, dp['number_of_chromosomes'] + 1):
      chr_shortname = 'chr' + str(c)
      lines_filepath = lines_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))
      genotype_filepath = genotype_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))
      variants_filepath = variants_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))

      locations.append(dict(cwd=args.working_directory, filetype='line', filename=lines_filepath))
      locations.append(dict(cwd=args.working_directory, filetype='genotype', filename=genotype_filepath))
      locations.append(dict(cwd=args.working_directory, filetype='variant', filename=variants_filepath))

    # Go through all the single files that are not named based off of a chromsome
    # Construct the file descriptor dictionaries, and then loop through and test each file's existance
    # phenotype_filename = Template('${cwd}/${growout}.ph.csv') # Welp, this is another instance of pheno file issue
    locations.append(dict(cwd=args.working_directory, filetype='kinship', filename=dp['kinship_filename']))
    locations.append(dict(cwd=args.working_directory, filetype='population_structure', filename=dp['population_structure_filename']))

    # Since there can be more than one file for the phenotypes, results, and run
    # For each array in the configuration file, add it to the list of paths to 
    # verify as existing

    for configuration_entry in dp:
      if isinstance(dp[configuration_entry], list):
        for filename in dp[configuration_entry]:
          locations.append(dict(cwd=args.working_directory, filetype=configuration_entry, filename=filename))
      else:
        # For any of the entries that CAN be a list, add their single values to
        # the file list
        if configuration_entry in [ 'phenotype_filename', 'gwas_run_filename', 'gwas_results_filename' ]:
          locations.append(dict(cwd=args.working_directory, filetype=configuration_entry, filename=dp[configuration_entry]))

    logging.debug("FILE LOCATIONS: %s", locations)

    for file_descriptor in locations:
      file_path = filepath_template.substitute(file_descriptor)
      if not os.path.isfile(file_path):
          raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


    logging.info(f'Found all files. Validating file contents.')

    # Validate the contents of each file
    for file_descriptor in locations:
      ft = file_descriptor['filetype']
      fp = filepath_template.substitute(file_descriptor)
      if ft == 'line':
        validate_line(conn, args, fp)
      elif ft == 'variant':
        validate_variant(conn, args, fp)
      elif ft == 'genotype':
        validate_genotype(conn, args, fp)
      elif ft == 'kinship':
        validate_kinship(conn, args, fp)
      elif ft == 'population_structure':
        validate_population_structure(conn, args, fp)
      elif ft == 'phenotype_filename':
        validate_phenotype(conn, args, fp)
      elif ft == 'gwas_run_filename':
        validate_runs(conn, args, fp)
      elif ft == 'gwas_results_filename':
        validate_results(conn, args, fp)
      else:
        logging.debug(f"Calling validation on unknown file: {fp}")
      
      logging.info(f"VALIDATED '{fp}'")
  except:
    raise

  logging.info(f'All input files appear to be valid.')


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
  parser.add_argument("-f", "--filename", action="store", help="Specify a configuration file. See documentation for expected format.")
  parser.add_argument("working_directory", action="store", metavar="WORKING_DIRECTORY", default=".", help="Working directory. Must contains all required files.")
  parser.add_argument("-c", "--config", default=default_config_location, help="JSON format file that contains metadata about the dataset to import.")
  parser.add_argument("--init", action="store_true", help="Scans directory for all files required to perform importation. Generates suggested JSON file used for importation")
  parser.add_argument("--design", action="store_true", help="Takes JSON file and validates if experiment can be imported")
  parser.add_argument("--collect", action="store_true", help="Gathers all the additional information in preperation for running GWAS on the data")
  parser.add_argument("--result", action="store_true", help="Collects last bit of user input to import the GWAS run and results information")
  parser.add_argument("--validate", action="store_true", help="Basic input file validate. Checks that all of the files specified in the configuration file exist and contain necessary data. Checks some basic data types.")
  parser.add_argument("--skip_genotype_validation", action="store_true", help="Errors in .012 files are infrequent, so enable this option to assume valid input.")

  args = parser.parse_args()
  if args.debug is True:
    args.verbose = True
    args.write = False

  if args.debug:
    pprint(args)
  
  logging_level = logging.INFO

  if args.debug:
    logging_level = logging.DEBUG
  
  logging_format = '%(asctime)s - %(levelname)s - %(filename)s %(lineno)d - %(message)s'
  logging.basicConfig(format=logging_format, level=logging_level)
    
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
    if args.validate:
      validate(args)
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
