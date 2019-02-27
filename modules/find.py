"""Looks up IDs of various elements in the database"""

import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
from dbconnect import config, connect
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result

def find_species(conn, speciesShortname):
  """Finds species by shortname 

    This function finds species_id for a species by its shortname

    :param conn: psycopg2 connection
    :type conn: connection object
    :param speciesShortname: human-readable shortname of species
    :type speciesShortname: string
    :return: species_id
    :rtype: integer
    
  """
  cur = conn.cursor()
  cur.execute("SELECT species_id FROM species WHERE shortname = '%s';" %  speciesShortname)
  row = cur.fetchone()
  if row is not None:
    speciesID = row[0]  
    cur.close()
    return speciesID
  else:
    return None

def find_population(conn, populationName):
  """Finds species by population name  

    This function finds the population_id for a population by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param populationName: human-readable name of population
    :type populationName: string
    :return: population_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT population_id FROM population WHERE population_name = '%s';" % populationName)
  row = cur.fetchone()
  if row is not None:
    populationID = row[0]
    cur.close()
    return populationID
  else:
    return None

def find_chromosome(conn, chromosome_name, chromosome_species):
  """Finds chromosome by name and species id

    This function finds the chromosome_id for a chromosome by its name and species id

    :param conn: psycopg2 connection
    :type conn: connection object
    :param chromosome_name: abbreviation of choromosome name
    :type chromosome_name: string
    :param chromosome_species: :ref:`species id <species_class>`
    :type chromosome_species: integer
    :return: chromosome_id
    :rtype: integer
  """
  cur = conn.cursor()
  # not sure if next line is correct...
  # TODO(timp): Check if this line meets functional requirements
  cur.execute("SELECT chromosome_id FROM chromosome WHERE chromosome_name = '%s' AND chromosome_species = '%s';" % (chromosome_name, chromosome_species))
  row = cur.fetchone()
  if row is not None:
    chromosomeID = row[0]
    cur.close()
    return chromosomeID
  else:
    return None

def find_line(conn, line_name, line_population):
  """Finds line by its name and population name 

    This function finds the line_id for a line by its name and population id

    :param conn: psycopg2 connection
    :type conn: connection object
    :param line_name: human-readable name of line
    :type line_name: string
    :param line_population: :ref:`population id <population_class>`
    :type line_population: integer
    :return: line_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT line_id FROM line WHERE line_name = '%s' AND line_population = '%s';" % (line_name, line_population))
  row = cur.fetchone()
  if row is not None:
    lineID = row[0]
    cur.close()
    return lineID
  else:
    return None

def find_growout_type(conn, growout_type):
  """Finds growout type by its name  

    This function finds the growout_id for a growout type by its type

    :param conn: psycopg2 connection
    :type conn: connection object
    :param growout_type: human-readable name of growout type
    :type growout_type: string
    :return: growout_type_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT growout_type_id FROM growout_type WHERE growout_type = '%s';" % growout_type)
  row = cur.fetchone()
  if row is not None:
    growout_type_ID = row[0]
    cur.close()
    return growout_type_ID
  else:
    return None

def find_growout(conn, growout_name):
  """Finds growout ID by its name
  
  This function finds the growout_id for a growout by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param growout_name: human-readable name of growout type
    :type growout_name: string
    :return: growout_id
    :rtype:integer
  """
  cur = conn.cursor()
  cur.execute(
      "SELECT growout_id FROM growout WHERE growout_name = '%s';" % growout_name)
  row = cur.fetchone()
  if row is not None:
    growout_id = row[0]
    cur.close()
    return growout_id
  else:
    return None
   
def find_location(conn, code):
  """Finds location by its code 

    This function finds the location_id for a location by its code

    :param conn: psycopg2 connection
    :type conn: connection object
    :param code: human-readable code assigned to a location
    :type code: string
    :return: location_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT location_id FROM location WHERE code = '%s';" % code)
  row = cur.fetchone()
  if row is not None:
    location_ID = row[0]
    cur.close()
    return location_ID
  else:
    return None

def find_kinship_algorithm(conn, algorithm):
  """Finds kinship algorithm by algorithm name 

    This function finds the kinship_algorithm_id by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param algorithm: human-readable name for algorithm
    :type algorithm: string
    :return: kinship_algorithm_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT kinship_algorithm_id FROM kinship_algorithm WHERE kinship_algorithm = '%s';" % algorithm)
  row = cur.fetchone()
  if row is not None:
    kinship_algorithm_ID = row[0]
    cur.close()
    return kinship_algorithm_ID
  else:
    return None

def find_population_structure_algorithm(conn, algorithm):
  """Finds population structure algorithm by algorithm name

    This function finds the population_structure_id by the name of its algorithm

    :param conn: psycopg2 connection
    :type conn: connection object
    :param algorithm: human-readable name for algorithm
    :type algorithm: string
    :return: population_structure_algorithm_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT population_structure_algorithm_id FROM population_structure_algorithm WHERE population_structure_algorithm = '%s';" % algorithm)
  row = cur.fetchone()
  if row is not None:
    population_structure_algorithm_ID = row[0]
    cur.close()
    return population_structure_algorithm_ID
  else:
    return None

def find_gwas_algorithm(conn, gwas_algorithm):
  """Finds algorithm used for genome-wide association study by algorithm name

    This function finds the gwas_algorithm_id by the name of the algorithm used in a genome-wide association study

    :param conn: psycopg2 connection
    :type conn: connection object
    :param gwas_algorithm: human-readable name of a GWAS algorithm
    :type gwas_algorithm: string
    :return: gwas_algorithm_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT gwas_algorithm_id FROM gwas_algorithm WHERE gwas_algorithm = '%s';" % gwas_algorithm)
  row = cur.fetchone()
  if row is not None:
    gwas_algorithm_ID = row[0]
    cur.close()
    return gwas_algorithm_ID
  else:
    return None

def find_genotype_version(conn, genotype_version_name):
  """Finds version of genotype by name 

    This function finds the genotype_version_id of a genotype by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param genotype_version_name: human-readable name of genotype version
    :type genotype_version_name: string
    :return: genotype_version_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT genotype_version_id FROM genotype_version WHERE genotype_version_name = '%s';" % genotype_version_name)
  row = cur.fetchone()
  if row is not None:
    genotype_version_ID = row[0]
    cur.close()
    return genotype_version_ID
  else:
    return None

def find_imputation_method(conn, imputation_method):
  """Finds imputation methodo by name

    This function finds the imputation_method_id by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param imputation_method: human-readable name of method
    :type imputation_method: string
    :return: imputation_method_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT imputation_method_id FROM imputation_method WHERE imputation_method = '%s';" % imputation_method)
  row = cur.fetchone()
  if row is not None:
    imputation_method_ID = row[0]
    cur.close()
    return imputation_method_ID
  else:
    return None

def find_kinship(conn, kinship_file_path):
  """Finds kinship by its location on a file system 

    This function finds the kinship_id by its location within a file system

    :param conn: psycopg2 connection
    :type conn: connection object
    :param kinship_file_path: path to kinship file
    :example path: ``/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv``
    :type kinship_file_path: string
    :return: kinship_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT kinship_id FROM kinship WHERE kinship_file_path = '%s';" % kinship_file_path)
  row = cur.fetchone()
  if row is not None:
    kinship_ID = row[0]
    cur.close()
    return kinship_ID
  else:
    return None

def find_population_structure(conn, population_structure_file_path):
  """Finds population_structure by its location within a file system 

    This function placeholder

    :param conn: psycopg2 connection
    :type conn: connection object
    :param population_structure_file_path: path to kinship file
    :example path: ``/opt/BaxDB/file_storage/population_structure_files/4.Eigenstrat.population.structure.10PCs.csv``
    :type population_structure_file_path: string
    :return: population_structure_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT population_structure_id FROM population_structure WHERE population_structure_file_path = '%s';" % population_structure_file_path)
  row = cur.fetchone()
  if row is not None:
    population_structure_ID = row[0]
    cur.close()
    return population_structure_ID
  else:
    return None

def find_trait(conn, trait_name):
  """Finds trait by its name 

    This function finds the traid_id for a trait by its name

    :param conn: psycopg2 connection
    :type conn: connection object
    :param traitName: human-readable name of trait
    :type traitName: string
    :return: trait_id
    :rtype: integer
  """
  cur = conn.cursor()
  cur.execute("SELECT trait_id FROM trait WHERE trait_name = '%s';" % trait_name)
  row = cur.fetchone()
  if row is not None:
    trait_ID = row[0]
    cur.close()
    return trait_ID
  else:
    return None

def find_gwas_run(conn, gwas_algorithm, missing_snp_cutoff_value, missing_line_cutoff_value, gwas_run_imputation_method, gwas_run_trait, nsnps, nlines, gwas_run_genotype_version, gwas_run_kinship, gwas_run_population_structure, minor_allele_frequency_cutoff_value):
  """Finds GWAS run by its parameters

    This function finds the gwas_run_id by its parameters

    :param conn: psycopg2 connection
    :type conn: connection object
    :param gwas_algorithm: :ref:`gwas_algorithm_id <gwas_algorithm_class>`
    :type gwas_algorithm: integer
    :param missing_snp_cutoff_value: ``todo``
    :type missing_snp_cutoff_value: numeric
    :param missing_line_cutoff_value: ``todo``
    :type missing_line_cutoff_value: numeric
    :param gwas_run_imputation_method: :ref:`imputation_method_id <imputation_method_class>`
    :type gwas_run_imputation_method: integer
    :param gwas_run_trait: :ref:`traid_id <trait_class>`
    :type gwas_run_trait: integer
    :param nsnps: ``todo``
    :type nsnps: integer
    :param nlines: ``todo``
    :type nlines: integer
    :param gwas_run_genotype_version: :ref:`genotype_version_id <genotype_version_class>`
    :type gwas_run_genotype_version: integer
    :param gwas_run_kinship: kinship id
    :type gwas_run_kinship: integer
    :param gwas_run_population_structure: :ref:`population_structure_id <population_structure_class>`
    :type gwas_run_population_structure: integer
    :param minor_allele_frequency_cutoff_value: ``todo``
    :type minor_allele_frequency_cutoff_value: numeric
    :return: gwas_run_id
    :rtype: integer

    .. note::
      Needs additional information on the
        - missing_snp_cutoff_value
        - missing_line_cutoff_value
        - nsnps
        - nlines
        - minor_allele_frequency_cutoff_value

  """
  cur = conn.cursor()
  cur.execute("SELECT gwas_run_id FROM gwas_run WHERE gwas_run_gwas_algorithm = '%s' AND missing_snp_cutoff_value = '%s' AND missing_line_cutoff_value = '%s' AND gwas_run_imputation_method = '%s' AND gwas_run_trait = '%s' AND nsnps = '%s' AND nlines = '%s' AND gwas_run_genotype_version = '%s' AND gwas_run_kinship = '%s' AND gwas_run_population_structure = '%s' AND minor_allele_frequency_cutoff_value = '%s';" % (gwas_algorithm, missing_snp_cutoff_value, missing_line_cutoff_value, gwas_run_imputation_method, gwas_run_trait, nsnps, nlines, gwas_run_genotype_version, gwas_run_kinship, gwas_run_population_structure, minor_allele_frequency_cutoff_value))
  row = cur.fetchone()
  if row is not None:
    gwas_run_ID = row[0]
    cur.close()
    return gwas_run_ID
  else:
    return None
