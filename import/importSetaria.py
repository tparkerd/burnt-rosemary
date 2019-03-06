#!/usr/bin/env python
"""Hard-coded importation script for Setaria Diverity Panel dataset.
The goal of this is to just get the data into the database and have a better grasp on how to generalize this

The MLMM GWAS Pipeline (by Greg Ziegler) is the source for the majority of information.

.. GitHub: https://github.com/gziegler/SetariaGWASPipeline

Curated Values & Files:
  The following values were confirmed by Greg to be accurate and correct to define the Setaria dataset.
  * Species shortname:
  * Species binomial name:
  * Species subspecies:
  * Species variety:
  * Population name:
  * Number of chromosome:
  * Lines filename (.indv): 
  * Genotype version name:
  * Genotype version: 
  * Reference genome line (line name):
  * Phenotype filename(s):


Todo:
  * Move the lines to a separate, *single* file because right now there are duplicates 
    of each said file for *each* chromosome.
  * Split the phenotypes/traits into individual files with `.ph` file extension

"""
import argparse
import csv
import sys
from string import Template
from random import randint
from pprint import pprint

import numpy as np
import pandas as pd
import psycopg2

from util import find, insert
from util.dbconnect import connect
from util.models import (chromosome, genotype, genotype_version, growout,
                         growout_type, gwas_algorithm, gwas_result, gwas_run,
                         imputation_method, kinship, kinship_algorithm, line,
                         location, phenotype, population, population_structure,
                         population_structure_algorithm, species, trait,
                         variant)
def process(args):
  """Workhorse function"""
  # =======================================
  # ========= Database Connection =========
  # =======================================
  try:
    conn = connect()
  except:
    raise

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
  species_shortname = 'sertaria' # setaria
  species_binomial = 'Setaria italica' # Setaria italica OR Setaria viridis  ???
  species_subspecies = None
  species_variety = None
  # Population
  population_name = 'SetariaPopulationName'
  # Chromosome
  chromosome_count = 9 # As defined by `awk -F'\t' '!a[$1]++{print NR":"$0}' 2.from12.setaria.maf0.1.maxMissing0.1.allLines.012.pos`
  # Line
  lines_filename = Template('../data/setaria/data/$chromosome_setaria.012.indv') # NOTE(tparker): Can use any chromosome, as they are the same for each. In the future, this the extraneous copies of the lines may be removed and there will be one specific line file, much like the phenotype files
  # Genotype Version
  # NOTE(tparker): This is possibly just the info about the reference genome
  #             It is likely included with the VCF genotype file (.012).
  genotype_version_name = 'SetariaGenotypeVersionName'
  genotype_version_identifier = 'SetariaGenotypeVersion' # NOTE(tparker): Not sure where to find this info or who names it
  line_name = 'SetariaReferenceGenomeLineName'
  # Growout, Type, and Location
  # NOTE(tparker): Unknown at this time
  ## Location
  ## Type
  ## Growout
  #
  # Traits
  phenotype_filename = Template('../data/setaria/$chromosome.ph.csv')

  # Model Construction & Insertion
  if not args.debug:
    # Species
    s = species(species_shortname, species_binomial, species_subspecies, species_variety)
    species_id = insert.insert_species(conn, s)
    # Population
    p = population(population_name, species_id)
    population_id = insert.insert_population(conn, p)
    # Chromosome
    chromosome_ids = insert.insert_all_chromosomes_for_species(conn, chromosome_count, species_id)
    # Line
    line_ids = insert.insert_lines_from_file(conn, lines_filename.safe_substitute('chr1'), population_id) # hard-coded substitue until just one file is used for lines
    # Genotype Version
    reference_genome_id = None # TODO(tparker): Look up line_id given reference genome (Not sure how this is selected), but it's found by its population_id (known) and a line_name
    line_id = find.find_line(conn, line_name, population_id) # NOTE(tparker): Need input from Greg on if this is the reference genome represented by a line
    gv = genotype_version(genotype_version_name,
                          genotype_version_identifier,
                          reference_genome = line_id,
                          genotype_version_population = population_id)
    genotype_version_id = insert.insert_genotype_version(conn, gv)

    # Growout, Type, and Location
    # NOTE(tparker): Unknown at this time
    ## Location
    ## Type
    ## Growout

    # Traits
    # Go through all the phenotype files available for the dataset and insert
    # the recorded traits for each.
    traits = list(pd.read_csv(phenotype_filename, index_col=0))
    trait_ids = insert.insert_traits_from_traitlist(conn, traits)

  # DEBUG
  else:
    print('Experiment Design\n=======================================')
    # Species
    s = species(species_shortname, species_binomial, species_subspecies, species_variety)
    print('\n------------------------\nSpecies\n------------------------')
    print(s)
    species_id = randint(1, 1000)
    print(f'Species ID set to {species_id}')
    # Population
    p = population(population_name, species_id)
    print('\n------------------------\nPopulation\n------------------------')
    print(p)
    population_id = randint(1, 1000)
    print(f'Population ID set to {population_id}')
    # Chromosome
    print('\n------------------------\nChromosome (from file)\n------------------------')
    print(f'insert_all_chromosomes_for_species(conn, {chromosome_count}, {species_id})')
    # Line
    print('\n------------------------\nLines (from file)\n------------------------')
    shortname_dict = dict(chromosome='chr1')
    print(f'insert_lines_from_file(conn, {lines_filename.safe_substitute(shortname_dict)}, {population_id})')
    # Genotype Version
    reference_genome_id = None
    line_id = randint(1, 1000)
    print('\n------------------------\nLine ID (Reference Genome)\n------------------------')
    print(f'Line ID set to {line_id}')
    gv = genotype_version(genotype_version_name,
                          genotype_version_identifier,
                          reference_genome = line_id,
                          genotype_version_population = population_id)
    print('\n------------------------\nGenotype Version\n------------------------')
    print(gv)
    genotype_version_id = randint(1, 1000)
    print(f'Genotype Version ID set to {genotype_version_id}')

    # Growout, Type, and Location
    # NOTE(tparker): Unknown at this time
    ## Location
    ## Type
    ## Growout

    # Traits
    # Go through all the phenotype files available for the dataset and insert
    # the recorded traits for each.
    
    print('\n------------------------\nTraits\n------------------------')
    print('Trait (from file)')

    print(f'list(pd.read_csv({phenotype_filename.safe_substitute(dict(chromosome="chr1"))}, index_col=0))')
    traits = [ 'weight', 'height', 'root_angle' ]
    print(f'insert.insert_traits_from_traitlist(conn, {traits})')
    trait_ids = [ randint(1, 1000) for t in traits ]
    print(f'Trait IDs set to {trait_ids}')

  # # =====================================
  # # ========== Pipeline Design ==========
  # # =====================================
  # # GWAS Algorithm: "MLMM", "EMMAx", "GAPIT", "FarmCPU"
  # # Imputation Method: "impute to major allele"
  # # Kinship Algorithm: "loiselle"
  # # Population Structure Algorithm: "Eigenstrat"

    # Expected User Input
    # GWAS Algorithm
    gwas_algorithm_name = 'MLMM' # According to Greg's README
    # Imputation Method
    imputation_method_name = 'SetariaImputationMethodName' # Unknown, apparently it was done by someone named Sujan
    # Kinship Algorithm
    kinship_algorithm_name = 'AstleBalding synbreed' # Placeholder, I don't know the exact string that should be used
    # Population Structure Algorithm
    population_structure_algorithm_name = 'Eigenstrat' # This is a guess based on filename

  if not args.debug:
    # Model Construction & Insertion
    # GWAS Algorithm
    gwas_algorithm_id = insert.insert_gwas_algorithm(conn, gwas_algorithm_name)
    # Imputation Method
    imputation_method_id = insert.insert_imputation_method(conn, imputation_method_name)
    # Kinship Algorithm
    kinship_algorithm_id = insert.insert_kinship_algorithm(conn, kinship_algorithm_name)
    # Population Structure Algorithm
    population_structure_algorithm_id = insert.insert_population_structure_algorithm(conn, population_structure_algorithm_name)

  else:
    print('\n\nPipeline Design\n=======================================')
    # Model Construction & Insertion
    # GWAS Algorithm
    print('\n------------------------\nGWAS Algorithm\n------------------------')
    print(f'insert_gwas_algorithm(conn, {gwas_algorithm_name})')
    gwas_algorithm_id = randint(1, 1000)
    print(f'GWAS Algorithm ID set to {gwas_algorithm_id}')
    # Imputation Method
    print('\n------------------------\Imputation Method\n------------------------')
    print(f'insert.insert_imputation_method(conn, {imputation_method_name})')
    imputation_method_id = randint(1, 1000)
    print(f'Imputation method ID set to {imputation_method_id}')
    # Kinship Algorithm
    print('\n------------------------\nKinship Algorithm\n------------------------')
    print(f'insert.insert_kinship_algorithm(conn, {kinship_algorithm_name})')
    kinship_algorithm_id = randint(1, 1000)
    print(f'Kinship algorithm ID set to {kinship_algorithm_id}')
    # Population Structure Algorithm
    print('\n------------------------\nPopulation Structure Algorithm\n------------------------')
    print(f'insert.insert_population_structure_algorithm(conn, {population_structure_algorithm_name})')
    population_structure_algorithm_id = randint(1, 1000)
    print(f'Population structure algorithm ID set to {population_structure_algorithm_id}')

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
  genotype_filename = Template('../data/setaria/$chromosome_shortname.012')
  # Variants
  variants_filename = Template('../data/setaria/$chromosome_shortname.012.pos')

  if not args.debug:
    # Model Construction & Insertion
    # Phenotype
    phenotype_ids = insert.insert_phenotypes_from_file(conn, phenotype_filename, population_id)
    # Genotype
    for c in range(1, chromosome_count + 1):
      chromosome_shortname = 'chr' + str(c)
      chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
      geno_filename = genotype_filename.substitute(chromosome_shortname)
      line_filename = lines_filename.substitute(chromosome_shortname)
      genotype_ids = insert.insert_genotypes_from_file(conn,
                                                      geno_filename,
                                                      line_filename,
                                                      chromosome_id,
                                                      population_id,
                                                      line_id)
    # Variants
    for c in range(1, chromosome_count + 1):
      chromosome_shortname = 'chr' + str(c)
      chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
      variant_filename = variants_filename.substitute(chromosome_shortname)
      variant_ids = insert.insert_variants_from_file(conn,
                                                    variant_filename,
                                                    species_id,
                                                    chromosome_id)

  else:
    print('\n\nExperiment Collection\n=======================================')
    # Model Construction & Insertion
    # Phenotype
    print('\n------------------------\nPhenotypes\n------------------------')
    print(f'insert.insert_phenotypes_from_file(conn, {phenotype_filename.safe_substitute(dict(chromosome="chr1"))}, {population_id})')
    # Genotype
    for c in range(1, chromosome_count + 1):
      chromosome_shortname = 'chr' + str(c)
      # chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
      chromosome_id = randint(1, 1000)
      print(f'Chromosome ID set to {chromosome_id}')
      geno_filename = genotype_filename.safe_substitute(dict(chromosome_shortname=chromosome_shortname))
      line_filename = lines_filename.safe_substitute(dict(chromosome_shortname=chromosome_shortname))
      print(f'insert.insert_genotypes_from_file(conn, {geno_filename}, {line_filename}, {chromosome_id}, {population_id}, {line_id})')
      genotype_ids = [ randint(1,1000) for g in range(1,25) ]
      print(f'Genotype IDs set to {genotype_ids}')
    # Variants
    for c in range(1, chromosome_count + 1):
      chromosome_shortname = 'chr' + str(c)
      # chromosome_id = find.find_chromosome(conn, chromosome_shortname, species_id)
      chromosome_id = randint(1, 1000)
      print(f'Chromosome ID set to {chromosome_id}')
      variant_filename = variants_filename.safe_substitute(dict(chromosome_shortname=chromosome_shortname))
      print(f'insert.insert_variants_from_file(conn, {variant_filename}, {species_id}, {chromosome_id})')


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
  kinship_filepath = '6.AstleBalding.synbreed.kinship.csv'
  # Population Structure
  # NOTE(tparker): Same reasoning as the kinship file. There should be a way 
  #                for the data to be stored in the database, not a 
  population_structure_filepath = '4.Eigenstrat.population.structure.10PCs.csv'

  if not args.debug:
    # Model Construction & Insertion
    # Kinship
    k = kinship(kinship_algorithm_id, kinship_filepath)
    kinship_id = insert.insert_kinship(conn, k)
    # Population Structure
    ps = population_structure(population_structure_algorithm_id, population_structure_filepath)
    population_structure_id = insert.insert_population_structure(conn, ps)
  else:
    print('\n\nPipeline Collection\n=======================================')
    # Model Construction & Insertion
    # Kinship
    k = kinship(kinship_algorithm_id, kinship_filepath)

    print('\n------------------------\nKinship\n------------------------')
    print(f'insert.insert_kinship(conn, {k})')
    kinship_id = randint(1, 1000)
    # Population Structure
    print('\n------------------------\nPopulation Structure Algorithm\n------------------------')
    ps = population_structure(population_structure_algorithm_id, population_structure_filepath)
    print(f'insert.insert_population_structure(conn, {ps})')
    population_structure_id = randint(1, 1000)
    print(f'Population structure ID set to {population_structure_id}')


  # =============================================
  # ================== Results ==================
  # =============================================
  # GWAS Run
  # GWAS Results

  # Expected User Input
  # GWAS Run
  gwas_run_filename = '../data/9.mlmmResults.csv'
  # The following values (0.2, 0.2, and 0.1) were all taken from the Maize282 import
  # NOTE(tparker): Make sure to double check with Greg on what the true values should be
  #                Also, double check the source of the pipeline to see if there is any
  #                indication what the values shoudl be.
  missing_snp_cutoff_value = 0.2
  missing_line_cutoff_value = 0.2
  minor_allele_frequency_cutoff_value = 0.1
  # GWAS Results
    # The following values (0.2, 0.2, and 0.1) were all taken from the Maize282 import
  # NOTE(tparker): Make sure to double check with Greg on what the true values should be
  #                Also, double check the source of the pipeline to see if there is any
  #                indication what the values shoudl be.
  gwas_results_filename = '../data/9.mlmmResults.csv'
  missing_snp_cutoff_value = 0.2
  missing_line_cutoff_value = 0.2
  minor_allele_frequency_cutoff_value = 0.1


  if not args.debug:
    # Model Construction & Insertion
    # GWAS Run
    # NOTE(tparker): Check with Greg on what the imputation method was used. I believe it was
    #                set by someone named Sujan because imputation was done beforehand
    imputation_method_id = find.find_imputation_method(conn, "impute to major allele")
    gwas_run_ids = insert.insert_gwas_runs_from_gwas_results_file(conn,
                                                                        gwas_run_filename,
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
                                                                gwas_results_filename,
                                                                gwas_algorithm_id,
                                                                missing_snp_cutoff_value,
                                                                missing_line_cutoff_value,
                                                                imputation_method_id,
                                                                reference_genome_id,
                                                                kinship_id,
                                                                population_structure_id,
                                                                minor_allele_frequency_cutoff_value)
  else:
    print('\n\nResults\n=======================================')
    # Model Construction & Insertion
    # GWAS Run
    # NOTE(tparker): Check with Greg on what the imputation method was used. I believe it was
    #                set by someone named Sujan because imputation was done beforehand

    print('\n------------------------\nImputation Method\n------------------------')
    # Imputation Method ID was already set in a previous set. If this was done at a
    # time, then it will have to be searched for in the database.
    print(f'Imputation method ID set to {imputation_method_id}')
    print('\n------------------------\nGWAS Run\n------------------------')
    print(f'insert.insert_gwas_runs_from_gwas_results_file(conn, {gwas_run_filename}, {gwas_algorithm_id}, {reference_genome_id}, {missing_snp_cutoff_value}, {missing_line_cutoff_value}, {minor_allele_frequency_cutoff_value}, {imputation_method_id}, {kinship_id}, {population_structure_id})')
    gwas_run_ids = [ randint(1, 1000) for g in range(1,15) ]
    print(f'GWAS run IDs set to {gwas_run_ids}')
    # GWAS Results
    print('\n------------------------\nGWAS Result\n------------------------')
    print(f'insert.insert_gwas_results_from_file(conn,{species_id},{gwas_results_filename},{gwas_algorithm_id},{missing_snp_cutoff_value},{missing_line_cutoff_value},{imputation_method_id},{reference_genome_id},{kinship_id},{population_structure_id},{minor_allele_frequency_cutoff_value})')
    gwas_result_ids = [ randint(1, 1000) for g in range(1,15) ]
    print(f'GWAS result IDs set to {gwas_result_ids}')

def parseOptions():
  """
  Function to parse user-provided options from terminal
  """
  description='Importation script specificaly for Setaria dataset.'
  # usage="python import.py -i ./data_file_directory [-o ./output_directory]"
  parser = argparse.ArgumentParser(description=description)
  parser.add_argument("--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("--debug", action="store_true", help="Enables --verbose and disables writes to disk")
  parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0-alpha')
  args = parser.parse_args()
  if args.debug is True:
    args.verbose = True
    args.write = False

  if args.debug:
    pprint(args)
  
  return args

if __name__ == "__main__":
  args = parseOptions()
  process(args)
