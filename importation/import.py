#!/usr/bin/env python
"""Hard-coded importation script for - dataset.
***This is effectively a template for hardcoding a dataset***

NOTICE: This file is currently being drafted to be a generic (whole) data import.
        This will *not* work for a partial dataset.

The goal of this is to just get the data into the database and have a better grasp on how
to generalize importation given a strict file structure.

The - GWAS Pipeline (by Greg Ziegler) is the source for the majority of information.

GitHub: -

Importation happens in five stages:
    Experiment Design          Pipeline Design
            |                         |
            |                         |
            V                         V
  Experiment Collection       Pipeline Collection
                   \             /
                    \           /
                     \         /
                      V       V
                       Results
  
  Although the experiment and pipeline pathways are not dependent on one another, they must
  both done before importing *any* GWAS results.

Required values & files categorized by stage:
  Experiment Design:
    * Species shortname
    * Species binomial name
    * Species subspecies
    * Species variety
    * Population name
    * Number of chromosomes
    * Lines filename (.indv) 
    * Genotype version assembly name
    * Genotype version annotation name 
    * Reference genome line name (line name)
    * Phenotype filename(s) (.ph.csv)*
  Pipeline Design:
    * GWAS algorithm name
    * Imputation method name
    * Kinship algortihm name
    * Population structure algorithm name
  Experiment Collection:
    * Genotype filename(s) (.012)
    * Variants filenames (.012.pos)
    * Phenotype filename(s)
  Pipeline Collection:
    * Kinship filename (.csv)
    * Population structure filename (.csv)
  Results:
    * GWAS run filename (.csv)
    * GWAS results filename (.csv)
    * Missing SNP cutoff value
    * Missing line cutoff value
    * Minor allele frequency cutoff value

Todo:
  * Move the lines to a separate, *single* file because right now there are duplicates 
    of each said file for *each* chromosome.
  * Split the phenotypes/traits into individual files with `.ph` file extension
  * Change genotype_version table so that fields reflect the change in the database
    This change will happen in the `importation.util.insert` module 

Additional Info:
  During the initial runs to import the data into a local database, these were the average
  execution times

  Species: < 1 sec
  Population: < 1 sec
  Lines: < 30 secs
  Genotype Version: < 1 sec
  Phenotypes: < 5 mins
  Genotypes: ~75 mins
  Variants: ~5.8 hours 

  Estimated total time: 

"""
import argparse
import configparser
import csv
import errno
import json
import os
import sys
from pprint import pprint
from random import randint
from string import Template

import numpy as np
import pandas as pd
import psycopg2

from importation.util import find, insert
from importation.util.dbconnect import connect
from importation.util.models import (chromosome, genotype, genotype_version,
                                     growout, growout_type, gwas_algorithm,
                                     gwas_result, gwas_run, imputation_method,
                                     kinship, kinship_algorithm, line,
                                     location, phenotype, population,
                                     population_structure,
                                     population_structure_algorithm, species,
                                     trait, variant)


def process(args):
  """Imports hardcoded values for Setaria database. Many items are placeholder values."""
  # =======================================
  # ========= Database Connection =========
  # =======================================
  try:
    conn = connect()
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
                        "minor_allele_frequency_cutoff_value"
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
                        "minor_allele_frequency_cutoff_value"
                      ]

    empty_fields = []
    for rf in required_fields:
      if not dp[rf]:
        empty_fields.append(rf)
    if empty_fields:
      raise KeyError(f'The following keys must be defined. Empty strings are not permitted. Please modify your JSON configuration: {empty_fields}')

    print(f'Configuration file is valid. Verifying that all files exist.')

    # Verify that all files exist
    # Lines
    lines_filename = Template('${cwd}/${chr}_${shortname}.012.indv')
    # Genotype
    genotype_filename = Template('${cwd}/${chr}_${shortname}.012')
    # Variants
    variants_filename = Template('${cwd}/${chr}_${shortname}.012.pos')

    for c in range(1, dp['number_of_chromosomes'] + 1):
      chr_shortname = 'chr' + str(c)
      lines_filepath = lines_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))
      genotype_filepath = genotype_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))
      variants_filepath = variants_filename.substitute(dict(cwd=args.working_directory, shortname=dp['species_shortname'], chr=chr_shortname))

      paths = [ lines_filepath, genotype_filepath, variants_filepath ]
      for file_path in paths:
        if not os.path.isfile(file_path):
          raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)

    # Go through all the single files that are not named based off of a chromsome
    # Construct the file descriptor dictionaries, and then loop through and test each file's existance
    # phenotype_filename = Template('${cwd}/${growout}.ph.csv') # Welp, this is another instance of pheno file issue
    generic_filename = Template('${cwd}/${filename}')
    locations = [
      dict(cwd=args.working_directory, filename=dp['kinship_filename']),
      dict(cwd=args.working_directory, filename=dp['population_structure_filename']),
      dict(cwd=args.working_directory, filename=dp['gwas_run_filename']),
      dict(cwd=args.working_directory, filename=dp['phenotype_filename'])
    ]

    for file_descriptor in locations:
      file_path = generic_filename.substitute(file_descriptor)
      if not os.path.isfile(file_path):
          raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)

  except:
    raise
  
  print(f'Found all files. Processing import.')

  # =======================================
  # ========== Experiment Design ==========
  # =======================================
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

  # Expected User Input
  # Species
  species_shortname = dp['species_shortname']
  species_binomial = dp['species_binomial_name']
  species_subspecies = dp['species_subspecies']
  species_variety = dp['species_variety']
  # Population
  population_name = dp['population_name']
  # Chromosome
  chromosome_count = dp['number_of_chromosomes']
  # Line
  # NOTE(tparker): Can use any chromosome, as they are the same for each.
  # In the future, this the extraneous copies of the lines may be removed
  # and there will be one specific line file, much like the phenotype files
  lines_filename = Template('${cwd}/${chr}_${shortname}.012.indv')
  # Genotype Version
  # NOTE(tparker): This is possibly just the info about the reference genome
  #                It is likely included with the VCF genotype file (.012).
  genotype_version_assembly_name = dp['genotype_version_assembly_name']
  genotype_version_annotation_name = dp['genotype_version_annotation_name']
  reference_genome_line_name = dp['reference_genome_line_name']
  # Growout, Type, and Location
  # NOTE(tparker): Unknown at this time
  ## Location
  ## Type
  ## Growout
  #
  # Traits
  phenotype_filename = f'{args.working_directory}/{dp["phenotype_filename"]}'

  # Model Construction & Insertion
  # Species
  s = species(species_shortname, species_binomial, species_subspecies, species_variety)
  species_id = insert.insert_species(conn, s)
  species_id = find.find_species(conn, species_shortname) # For idempotence
  # Population
  p = population(population_name, species_id)
  population_id = insert.insert_population(conn, p)
  population_id = find.find_population(conn, population_name) # For idempotence
  if args.verbose:
    print(f'[Insert]\tPopulation ID\t{population_id}')
  # Chromosome
  chromosome_ids = insert.insert_all_chromosomes_for_species(conn, chromosome_count, species_id)
  # Line
  working_filepath = lines_filename.substitute(dict(chr="chr1", cwd=f"{args.working_directory}", shortname=species_shortname))
  try:
    if not os.path.isfile(working_filepath):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), working_filepath)
  except:
    raise

  line_ids = insert.insert_lines_from_file(conn, working_filepath, population_id) # hard-coded substitue until just one file is used for lines
  # Genotype Version
  reference_genome_id = find.find_line(conn, reference_genome_line_name, population_id)
  gv = genotype_version(genotype_version_assembly_name,
                        genotype_version_annotation_name,
                        reference_genome = reference_genome_id,
                        genotype_version_population = population_id)
  genotype_version_id = insert.insert_genotype_version(conn, gv)
  genotype_version_id = find.find_genotype_version(conn, genotype_version_assembly_name)
  
  # Growout, Type, and Location
  # NOTE(tparker): Unknown at this time
  ## Location
  ## Type
  ## Growout

  # Traits
  # Go through all the phenotype files available for the dataset and insert
  # the recorded traits for each.
  try:
    if not os.path.isfile(phenotype_filename):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), phenotype_filename)
  except:
    raise
  traits = list(pd.read_csv(phenotype_filename, index_col=0))
  trait_ids = insert.insert_traits_from_traitlist(conn, traits)

  # # =====================================
  # # ========== Pipeline Design ==========
  # # =====================================
  # # GWAS Algorithm: "MLMM", "EMMAx", "GAPIT", "FarmCPU"
  # # Imputation Method: "impute to major allele"
  # # Kinship Algorithm: "loiselle"
  # # Population Structure Algorithm: "Eigenstrat"

  # Expected User Input
  # GWAS Algorithm
  gwas_algorithm_name = dp['gwas_algorithm_name'] # According to Greg's README
  # Imputation Method
  imputation_method_name = dp['imputation_method_name'] # Unknown, apparently it was done by someone named Sujan
  # Kinship Algorithm
  kinship_algorithm_name = dp['kinship_algortihm_name'] # Placeholder, I don't know the exact string that should be used
  # Population Structure Algorithm
  population_structure_algorithm_name = dp['population_structure_algorithm_name'] # This is a guess based on filename

  # Model Construction & Insertion
  # GWAS Algorithm
  ga = gwas_algorithm(gwas_algorithm_name)
  gwas_algorithm_id = insert.insert_gwas_algorithm(conn, ga)
  gwas_algorithm_id = find.find_gwas_algorithm(conn, 'MLMM')
  # Imputation Method
  im = imputation_method(imputation_method_name)
  imputation_method_id = insert.insert_imputation_method(conn, im)
  # Kinship Algorithm
  ka = kinship_algorithm(kinship_algorithm_name)
  kinship_algorithm_id = insert.insert_kinship_algorithm(conn, ka)
  # Population Structure Algorithm
  psa = population_structure_algorithm(population_structure_algorithm_name)
  population_structure_algorithm_id = insert.insert_population_structure_algorithm(conn, psa)

  # ===========================================
  # ========== Experiment Collection ==========
  # ===========================================
  # Phenotype (external source?)
  #       This needs to be standardized to a .pheno filetype.
  #       For now, it is the longFormat for the Maize282 datset
  #       5.mergedWeightNorm.LM.rankAvg.longFormat.csv, but for Setaria will be 
  # Genotype (VCF output)
  # Variant (VCF output)

  # Expected User Input
  # Phenotype
  # NOTE(tparker): Define in earlier stage
  # Genotype
  genotype_filename = Template('${cwd}/${chr}_${shortname}.012')
  # Variants
  variants_filename = Template('${cwd}/${chr}_${shortname}.012.pos')

  # Model Construction & Insertion
  # Phenotype
  try:
    if not os.path.isfile(phenotype_filename):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), phenotype_filename)
  except:
    raise
  phenotype_ids = insert.insert_phenotypes_from_file(conn, phenotype_filename, population_id)
  # Genotype
  for c in range(1, chromosome_count + 1):
    chromosome_shortname = 'chr' + str(c)
    chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
    geno_filename = genotype_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.working_directory}', shortname=species_shortname))
    line_filename = lines_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.working_directory}', shortname=species_shortname))
    try:
      if not os.path.isfile(geno_filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), geno_filename)
      if not os.path.isfile(line_filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), line_filename)
    except:
      raise
    genotype_ids = insert.insert_genotypes_from_file(conn, geno_filename, line_filename, chromosome_id, population_id, genotype_version_id)
  # Variants
  for c in range(1, chromosome_count + 1):
    chromosome_shortname = 'chr' + str(c)
    chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
    variant_filename = variants_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.working_directory}', shortname=species_shortname))
    try:
      if not os.path.isfile(variant_filename):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), variant_filename)
    except:
      raise
    variant_ids = insert.insert_variants_from_file(conn,
                                                  variant_filename,
                                                  species_id,
                                                  chromosome_id)

  # =========================================
  # ========== Pipeline Collection ==========
  # =========================================
  # Kinship
  # Setaria Kinship is stored in:
  ## /shares/ibaxter_share/gziegler/SetariaGWASPipeline/data/genotype/6.AstleBalding.synbreed.kinship.rda 
  ## Exported the file to CSV using R
  ### load('6.AstleBalding.synbreed.kinship.rda')
  ### write.csv(kinship, '6.AstleBalding.synbreed.kinship.csv')
  # Population Structure

  # Expected User Input
  # Kinship
  # NOTE(tparker): Currently the database just stores the filename.
  #                There is no service to upload the file to database's
  #                host, so there's no single location to find the file
  #                I would like to find out why this is the case and if 
  #                it would just be better to store it in the database and
  #                allow the user to export the table themselves as a CSV.
  kinship_filepath = f'{args.working_directory}/{dp["kinship_filename"]}'
  # Population Structure
  # NOTE(tparker): Same reasoning as the kinship file. There should be a way 
  #                for the data to be stored in the database, not a 
  population_structure_filepath = f'{args.working_directory}/{dp["population_structure_filename"]}'

  # Model Construction & Insertion
  # Kinship
  try:
    if not os.path.isfile(kinship_filepath):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), kinship_filepath)
    if not os.path.isfile(population_structure_filepath):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), population_structure_filepath)
  except:
    raise
  k = kinship(kinship_algorithm_id, kinship_filepath)
  kinship_id = insert.insert_kinship(conn, k)
  # Population Structure
  ps = population_structure(population_structure_algorithm_id, population_structure_filepath)
  population_structure_id = insert.insert_population_structure(conn, ps)
  population_structure_id = find.find_population_structure(conn, population_structure_filepath)


  # =============================================
  # ================== Results ==================
  # =============================================
  # GWAS Run
  # GWAS Results

  # Expected User Input
  # GWAS Run & results
  gwas_filename = f'{args.working_directory}/{dp["gwas_results_filename"]}'
  # The following values (0.2, 0.2, and 0.1) were all taken from the Maize282 import
  # NOTE(tparker): Make sure to double check with Greg on what the true values should be
  #                Also, double check the source of the pipeline to see if there is any
  #                indication what the values shoudl be.
  missing_snp_cutoff_value = dp['missing_SNP_cutoff_value']
  missing_line_cutoff_value = dp['missing_line_cutoff_value']
  minor_allele_frequency_cutoff_value = dp['minor_allele_frequency_cutoff_value']

  # Model Construction & Insertion
  # GWAS Run
  # NOTE(tparker): Check with Greg on what the imputation method was used. I believe it was
  #                set by someone named Sujan because imputation was done beforehand
  try:
    if not os.path.isfile(gwas_filename):
      raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), gwas_filename)
  except:
    raise
  imputation_method_id = find.find_imputation_method(conn, imputation_method_name)
  gwas_run_ids = insert.insert_gwas_runs_from_gwas_results_file(conn,
                                                                gwas_filename,
                                                                gwas_algorithm_id,
                                                                reference_genome_id,
                                                                missing_snp_cutoff_value,
                                                                missing_line_cutoff_value,
                                                                minor_allele_frequency_cutoff_value,
                                                                imputation_method_id,
                                                                kinship_id,
                                                                population_structure_id)
  # GWAS Results
  gwas_result_ids = insert.insert_gwas_results_from_file(conn,
                                                         species_id,
                                                         gwas_filename,
                                                         gwas_algorithm_id,
                                                         missing_snp_cutoff_value,
                                                         missing_line_cutoff_value,
                                                         imputation_method_id,
                                                         reference_genome_id,
                                                         kinship_id,
                                                         population_structure_id,
                                                         minor_allele_frequency_cutoff_value)

def parseOptions():
  """
  Function to parse user-provided options from terminal
  """
  description="""Importation script specificaly for Setaria dataset.
  
  Importation happens in five stages:
    Experiment Design          Pipeline Design
            |                         |
            |                         |
            V                         V
  Experiment Collection       Pipeline Collection
                   \             /
                    \           /
                     \         /
                      V       V
                       Results

  Example configuration file:
    # configuration.json
    {
      "species": {
        "shortname": "setaria",
        "binomial name": "Setaria Dkss"
      },
      "population": "wg393",
      "growout": [
        {
          "name": "PH18",
          "year": "2018",
          "type": "phenotyper",
          "location": {
            "code": "PH18"
          }
        }
      ]
    }
    
    """
  parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("--debug", action="store_true", help="Enables --verbose and disables writes to disk")
  parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0-alpha')
  parser.add_argument("-f", "--filename", action="store", help="Specify a configuration file. See documentation for expected format.")
  parser.add_argument("--log", action="store_true", help="Enabled logging. Filename is appended to %(prog)s.log")
  parser.add_argument("working_directory", action="store", metavar="WORKING_DIRECTORY", default=".", help="Working directory. Must contains all required files.")
  args = parser.parse_args()
  if args.debug is True:
    args.verbose = True
    args.write = False

  if args.log is True:
    args.log_file = f'{os.path.basename(__file__)}'

  if args.debug:
    pprint(args)
  
  return args

if __name__ == "__main__":
  args = parseOptions()
  process(args)
