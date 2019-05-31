"""Constructs and executes SQL to insert data into database"""
import time

import numpy as np
import pandas as pd
import psycopg2 as pg
from tqdm import tqdm

import util.find as find
import util.parsinghelpers as ph
from util.models import (chromosome, genotype, genotype_version, growout,
                         growout_type, gwas_algorithm, gwas_result, gwas_run,
                         imputation_method, kinship, kinship_algorithm, line,
                         location, phenotype, population, population_structure,
                         population_structure_algorithm, species, trait,
                         variant)


def exists_in_database(cur, SQL):
  """Checks if an object has already been inserted into the database

  Args:
    cur (psycopg2.extensions.cursor): psycopg2 cursor
    SQL (str): SQL statement to be executed

  Returns:
    int if true, None otherwise.

  """
  cur.execute(SQL)
  query_result = cur.fetchone()
  if query_result is not None:
    known_id = query_result[0]
    if known_id:
      return known_id
  
  return None

def insert_species(conn, args, species):
  """Inserts species into database by its shortname, binomial, subspecies, and variety

  This function inserts a species into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    species (species): :ref:`species <species_class>` object
  
  Returns:
    int: species id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT species_id \
          FROM species \
          WHERE shortname = '{species.n}' AND \
                binomial = '{species.b}' AND \
                subspecies = '{species.s}' AND \
                variety = '{species.v}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO species (shortname, binomial, subspecies, variety)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING species_id;"""
  args_tuple = (species.n, species.b, species.s, species.v)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_population(conn, args, population):
  """Inserts population into database

  This function inserts a population into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    population (population): :ref:`population <population_class>` object

  Returns:
    int: population id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT population_id \
          FROM population \
          WHERE population_name = '{population.n}' AND \
                population_species = {population.s}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO population (population_name, population_species)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING population_id;"""
  args_tuple = (population.n, population.s)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
     newID = row[0]
     conn.commit()
     cur.close()
     return newID
  else:
    return None  

def insert_chromosome(conn, args, chromosome):
  """Inserts chromosome into database by its name

  This function inserts a chromosome into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    chromosome (chromosome): :ref:`chromosome <chromosome_class>` object
  
  Returns:
    int: chromosome id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT chromosome_id \
          FROM chromosome \
          WHERE chromosome_name = '{chromosome.n}' AND \
                chromosome_species = {chromosome.s}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    conn.commit()
    cur.close()
    return known_id


  SQL = """INSERT INTO chromosome (chromosome_name, chromosome_species)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING chromosome_id;"""
  args_tuple = (chromosome.n, chromosome.s)
  try:
    cur.execute(SQL, args_tuple)
  except pg.Error as err:
    print("%s: %s" % (err.__class__.__name__, err))
    raise
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_all_chromosomes_for_species(conn, args, numChromosomes, speciesID):
  """Inserts all chromosomes for a species into database by its name

  This function inserts all chromosomes for a species into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    numChromosomes (int): upper-bound number of chromosomes to consider for a species
    species (species): :ref:`species <species_class>` object
  
  Returns:
    list of int: list of species ids
  """
  chrlist = ph.generate_chromosome_list(numChromosomes)
  insertedChromosomeIDs = []
  for chrname in chrlist:
    chrobj = chromosome(chrname, speciesID)
    insertedChromosomeID = insert_chromosome(conn, args, chrobj)
    insertedChromosomeIDs.append(insertedChromosomeID)
  return insertedChromosomeIDs


def insert_line(conn, args, line):
  """Inserts line into database

  This function inserts a line into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    line (line): :ref:`line <line_class>` object
    
  Returns:
    int: line id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT line_id \
          FROM line \
          WHERE line_name = '{line.n}' AND \
                line_population = {line.p}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO line (line_name, line_population)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING line_id;"""
  args_tuple = (line.n, line.p)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_lines_from_file(conn, args, lineFile, populationID):
  """Inserts lines into database from a file

  This function inserts a lines into a database from a file


  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    lineFile (str): absolute path to input file
    populationID (int): :ref:`population <population_class>`
  
  Returns:
    list of int: list of population id
  """
  linelist = ph.parse_lines_from_file(lineFile)
  insertedLineIDs = []
  for linename in tqdm(linelist, desc="Lines"):
    lineobj = line(linename, populationID)
    insertedLineID = insert_line(conn, args, lineobj)
    insertedLineIDs.append(insertedLineID)
  return insertedLineIDs


def insert_variant(conn, args, variant):
  """Inserts variant into database

  This function inserts a variant into a database

  conn (psycopg2.extensions.connection): psycopg2 connection
  args (ArgumentParser namespace): user-defined arguments
  
  Args:
    variant (variant): :ref:`variant <variant_class>` object
  
  Returns:
    int: variant id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT variant_id \
          FROM variant \
          WHERE variant_species = {variant.s} AND \
                variant_chromosome = {variant.c} AND \
                variant_pos = {variant.p}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Variant found] {variant}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO variant(variant_species, variant_chromosome, variant_pos)
        VALUES (%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING variant_id;"""
  args_tuple = (variant.s, variant.c, variant.p)
  cur.execute(SQL, args_tuple)
  #newID = cur.fetchone()[0]
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_variants_from_file(conn, args, variantPosFile, speciesID, chromosomeID):
  """Inserts chromosome into database by its name

  This function inserts a chromosome into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    variantPosFile (str): absolute path to input file
    speciesID (int): :ref:`species <species_class>`
    chromosomeID (int): :ref:`chromosome <chromosome_class>`

  Returns:
    list of int: list of variant id
  """
  variantlist = ph.parse_variants_from_file(variantPosFile)
  cVariants = len(variantlist)
  insertedVariantIDs = []
  for variantpos in tqdm(variantlist, desc="Variants from %s" % variantPosFile):
    variantobj = variant(speciesID, chromosomeID, variantpos)
    insertedVariantID = insert_variant(conn, args, variantobj)
    insertedVariantIDs.append(insertedVariantID)
  return insertedVariantIDs


def insert_genotype(conn, args, genotype):
  """Inserts genotype into database

  This function inserts a genotype into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    genotype (genotype): :ref:`genotype <genotype_class>` object

  Returns:
    int: genotype id
  """
  cur = conn.cursor()

  # See if the genotype has already been inserted, and if so, return it
  SQL = f"SELECT genotype_id \
          FROM genotype \
          WHERE genotype_line = {genotype.l} AND \
                genotype_chromosome = {genotype.c} AND \
                genotype_genotype_version = '{genotype.v}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Genotype found] ({genotype.l}, {genotype.c}, {genotype.v})')
    conn.commit()
    cur.close()
    return known_id

  # Otherwise, the genotype has not been inserted yet, so do so
  SQL = """INSERT INTO genotype(genotype_line, genotype_chromosome, genotype, genotype_genotype_version)
        VALUES (%s,%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING genotype_id;"""

  args_tuple = (genotype.l, genotype.c, genotype.g, genotype.v)
  cur.execute(SQL, args_tuple)
  genotype_id = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return genotype_id


def insert_genotypes_from_file(conn, args, genotypeFile, lineFile, chromosomeID, populationID, genotype_versionID):
  """Inserts genotypes into database

  This function inserts a genotypes into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    genotypeFile (str): absolute path to input file
    lineFile (str): absolute path to input file
    chromosomeID (int): :ref:`chromosome <chromosome_class>`
    populationID (int): :ref:`population <population_class>`
  
  Returns:
    list of int: list of genotype id
  """
  genotypes = ph.parse_genotypes_from_file(genotypeFile)
  linelist = ph.parse_lines_from_file(lineFile)
  lineIDlist = ph.convert_linelist_to_lineIDlist(conn, args, linelist, populationID)
  zipped = zip(lineIDlist, genotypes)
  ziplist = list(zipped)
  insertedGenotypeIDs = []
  for zippedpair in tqdm(ziplist, desc="Genotypes from %s" % genotypeFile):
    genotypeObj = genotype(zippedpair[0], chromosomeID, zippedpair[1], genotype_versionID)
    insertedGenotypeID = insert_genotype(conn, args, genotypeObj)
    insertedGenotypeIDs.append(insertedGenotypeID)


  return insertedGenotypeIDs

def insert_growout(conn, args, growout):
  """Inserts growout into database

  This function inserts a growout into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    growout (growout): :ref:`growout <genotype_class>` object
    
  Returns:
    int: growout id
  """
  cur = conn.cursor()
  SQL = """INSERT INTO growout(growout_name, growout_population, growout_location, year, growout_growout_type)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING growout_id;"""
  args_tuple = (growout.n, growout.p, growout.l, growout.y, growout.t)
  try:
    cur.execute(SQL, args_tuple)
  except pg.Error as err:
    print("%s: %s" % (err.__class__.__name__, err))
    raise
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_location(conn, args, location):
  """Inserts location into database

  This function inserts a location into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    location (location): :ref:`location <location_class>` object
  
  Returns:
    int: location id
  """
  cur = conn.cursor()
  SQL = """INSERT INTO location(country, state, city, code)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING location_id;"""
  args_tuple = (location.c, location.s, location.i, location.o)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None  

def insert_phenotype(conn, args, phenotype):
  """Inserts phenotype into database

  This function inserts a phenotype into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    phenotype (phenotype): :ref:`phenotype <phenotype_class>` object
    
  Returns:
    int: phenotype id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT phenotype_id \
          FROM phenotype \
          WHERE phenotype_line = {phenotype.l} AND \
                phenotype_trait = {phenotype.t} AND \
                LOWER(phenotype_value) = '{phenotype.v}'" # Lower needed to deal with SQL's representation of NaN ('NaN') and Python's ('nan')
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Phenotype found] {phenotype}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO phenotype(phenotype_line, phenotype_trait, phenotype_value)
        VALUES (%s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING phenotype_id;"""
  args_tuple = (phenotype.l, phenotype.t, phenotype.v)
  try:
    cur.execute(SQL, args_tuple)
  except pg.Error as err:
    print("%s: %s" % (err.__class__.__name__, err))
    raise
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_phenotypes_from_file(conn, args, phenotype_filename, population_id):
  """Inserts phenotypes into database

  This function inserts phenotypes from a file into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    phenotype_filename (str): absolute path to input file
    population_id (int): :ref:`population_id <population_class>`

  Returns:
    list of int: list of phenotype id
  """
  # Read through just the first column of the CSV, ignoring any column
  df = pd.read_csv(phenotype_filename, index_col=0)
  phenotype_ids = []
  for key, value in tqdm(df.iteritems(), total=len(df.columns), desc="Phenotypes"):
    trait_id = find.find_trait(conn, args, key)
    for line_name, traitval in value.iteritems():
      line_id = find.find_line(conn, args, line_name, population_id)
      if line_id is None:
        l = line(line_name, population_id)
        line_id = insert_line(conn, args, l)
      ph = phenotype(line_id, trait_id, traitval)
      phenotype_id = insert_phenotype(conn, args, ph)
      phenotype_ids.append(phenotype_id)
  return phenotype_ids


def insert_trait(conn, args, trait):
  """Inserts trait into database

  This function inserts a trait into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    trait (trait): :ref:`trait <trait_class>` object

  Returns:
    int: trait id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT trait_id \
          FROM trait \
          WHERE trait_name = '{trait.n}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Trait found] {trait}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO trait(trait_name)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING trait_id;"""
  arg = (trait.n,)
  cur.execute(SQL, arg)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_traits_from_traitlist(conn, args, traitlist):
  """Inserts traits from list into database

  This function inserts a traitlist into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    traitlist (list of str): list of trait names
  
  Returns:
    list of int: list of trait id
  """
  traitIDs = []
  for traitname in tqdm(traitlist, desc="Traits"):
    traitObj = trait(traitname, None, None, None)
    insertedTraitID = insert_trait(conn, args, traitObj)
    traitIDs.append(insertedTraitID)
  return traitIDs


def insert_growout_type(conn, args, growout_type):
  """Inserts growout type into database

  This function inserts a growout type into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    growout_type (growout_type): :ref:`growout_type <growout_type_class>` object

  Returns:
    int: growout type id
  """
  cur = conn.cursor()
  SQL = """INSERT INTO growout_type(growout_type)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING growout_type_id;"""
  arg = (growout_type.t,)
  cur.execute(SQL, arg)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_gwas_algorithm(conn, args, gwas_algorithm):
  """Inserts GWAS algorithm into database

  This function inserts a GWAS algorithm into a database


  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    gwas_algorithm (gwas_algorithm): :ref:`gwas_algorithm <gwas_algorithm_class>` object
  
  Returns:
    int: gwas algorithm id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT gwas_algorithm_id \
          FROM gwas_algorithm \
          WHERE gwas_algorithm = '{gwas_algorithm.a}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Genotype algorithm found] {gwas_algorithm}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO gwas_algorithm(gwas_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_algorithm_id;"""
  args = (gwas_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_genotype_version(conn, args, genotype_version):
  """Inserts genotype version into database

  This function inserts a genotype version into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    genotype_version (genotype_version): :ref:`genotype_version <genotype_version_class>` object

  Returns:
    int: genotype version id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT genotype_version_id \
          FROM genotype_version \
          WHERE genotype_version_assembly_name = '{genotype_version.n}' AND \
                genotype_version_annotation_name = '{genotype_version.v}' AND \
                reference_genome = {genotype_version.r} AND \
                genotype_version_population = {genotype_version.p}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Genotype version found] {genotype_version}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO genotype_version(genotype_version_assembly_name, genotype_version_annotation_name, reference_genome, genotype_version_population)
        VALUES (%s,%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING genotype_version_id;"""
  # print("Genotype Version: " + str(genotype_version))
  args_tuple = (genotype_version.n, genotype_version.v, genotype_version.r, genotype_version.p)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_imputation_method(conn, args, imputation_method):
  """Inserts imputation method into database

  This function inserts a imputation method into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    imputation_method (imputation_method): :ref:`imputation_method <imputation_method_class>` object
  
  Returns:
    int: imputation method id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT imputation_method_id \
          FROM imputation_method \
          WHERE imputation_method = '{imputation_method.m}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Imputation method found] {imputation_method}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO imputation_method(imputation_method)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING imputation_method_id;"""
  args = (imputation_method.m,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_kinship_algorithm(conn, args, kinship_algorithm):
  """Inserts kinship_algorithm into database

  This function inserts a kinship_algorithm into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    kinship_algorithm (kinship_algorithm): :ref:`kinship_algorithm <kinship_algorithm_class>` object

  Returns:
    int: kinship algorithm id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT kinship_algorithm_id \
          FROM kinship_algorithm \
          WHERE kinship_algorithm = '{kinship_algorithm.a}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Kinship algorithm found] {kinship_algorithm}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO kinship_algorithm(kinship_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING kinship_algorithm_id;"""
  args = (kinship_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_kinship(conn, args, kinship):
  """Inserts kinship into database

  This function inserts a kinship into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    kinship (kinship): :ref:`kinship <kinship_class>` object
    
  Returns:
    int: kinship id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT kinship_id \
          FROM kinship \
          WHERE kinship_algorithm = '{kinship.a}' AND \
                kinship_file_path = '{kinship.p}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Kinship found] {kinship}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO kinship(kinship_algorithm, kinship_file_path)
        VALUES (%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING kinship_id;"""
  args = (kinship.a, kinship.p)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_population_structure_algorithm(conn, args, population_structure_algorithm):
  """Inserts population_structure_algorithm into database

  This function inserts a population_structure_algorithm into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    population_structure_algorithm (population_structure_algorithm): :ref:`population_structure_algorithm <population_structure_algorithm_class>` object
  
  Returns:
    int: population structure algorithm id
  """
  cur = conn.cursor()
  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT population_structure_algorithm_id \
          FROM population_structure_algorithm \
          WHERE population_structure_algorithm = '{population_structure_algorithm.a}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Population structure algorithm found] {population_structure_algorithm}')
    conn.commit()
    cur.close()
    return known_id
  
  SQL = """INSERT INTO population_structure_algorithm(population_structure_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING population_structure_algorithm_id;"""
  args = (population_structure_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_population_structure(conn, args, population_structure):
  """Inserts population into database

  This function inserts a population into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    population (population): :ref:`population <population_class>` object
  
  Returns:
    int: population id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT population_structure_id \
          FROM population_structure \
          WHERE population_structure_algorithm = '{population_structure.a}' AND \
                population_structure_file_path = '{population_structure.p}'"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[Population structure found] {population_structure}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO population_structure(population_structure_algorithm, population_structure_file_path)
        VALUES (%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING population_structure_id;"""
  args = (population_structure.a, population_structure.p)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_gwas_run(conn, args, gwas_run):
  """Inserts gwas_run into database

  This function inserts a gwas_run into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    gwas_run (gwas_run): :ref:`gwas_run <gwas_run_class>` object
  
  Returns:
    int: gwas run id
  """
  cur = conn.cursor()

  # See if data has already been inserted, and if so, return it
  SQL = f"SELECT gwas_run_id \
          FROM gwas_run \
          WHERE gwas_run_trait = {gwas_run.t} AND \
                nsnps = {gwas_run.s} AND \
                nlines = {gwas_run.l} AND \
                gwas_run_gwas_algorithm = {gwas_run.a} AND \
                gwas_run_genotype_version = {gwas_run.v} AND \
                missing_snp_cutoff_value = {gwas_run.m} AND \
                missing_line_cutoff_value = {gwas_run.i} AND \
                minor_allele_frequency_cutoff_value = {gwas_run.n} AND \
                gwas_run_imputation_method = {gwas_run.p} AND \
                gwas_run_kinship = {gwas_run.k} AND \
                gwas_run_population_structure = {gwas_run.o}"
  
  known_id = exists_in_database(cur, SQL)
  if known_id is not None:
    if args.verbose is True:
      print(f'[GWAS run found] {gwas_run}')
    conn.commit()
    cur.close()
    return known_id

  SQL = """INSERT INTO gwas_run(gwas_run_trait, nsnps, nlines, gwas_run_gwas_algorithm, gwas_run_genotype_version, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwas_run_imputation_method, gwas_run_kinship, gwas_run_population_structure)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_run_id;"""
  args = (gwas_run.t, gwas_run.s, gwas_run.l, gwas_run.a, gwas_run.v, gwas_run.m, gwas_run.i, gwas_run.n, gwas_run.p, gwas_run.k, gwas_run.o)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_gwas_runs_from_gwas_results_file(conn, args, gwas_results_file, gwasRunAlgorithmID, gwasRunGenotypeVersionID, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwasRunImputationMethodID, gwasRunKinshipID, gwasRunPopulationStructureID):
  """Inserts a collection of GWAS runs from an input file into database

  This function inserts a a collection of GWAS runs from an input file into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    gwas_results_file (str): absolute path to input file
    gwasRunAlgorithmID (int): :ref:`gwas_algorithm_id <gwas_algorithm_class>`
    gwasRunGenotypeVersionID (int): :ref:`genotype_version_id <genotype_version_class>`
    missing_snp_cutoff_value (numeric): 
    missing_line_cutoff_value (numeric):
    minor_allele_frequency_cutoff_value (numeric):
    gwasRunImputationMethodID (int): :ref:`imputation_method_id <imputation_method_class>`
    gwasRunKinshipID (int): :ref:`kinship_id <kinship_class>`
    gwasRunPopulationStructureID (int): :ref:`population_structure_id <population_structure_class>`
  Returns:
    list of int: GWAS Run IDs
  """
  gwas_run_list = ph.parse_unique_runs_from_gwas_results_file(gwas_results_file)
  insertedGwasRunIDs = []
  for gwas_run_item in tqdm(gwas_run_list, desc="GWAS Runs"):
    traitID = find.find_trait(conn, args, gwas_run_item[0])
    gwas_run_obj = gwas_run(traitID, gwas_run_item[1], gwas_run_item[2], gwasRunAlgorithmID, gwasRunGenotypeVersionID, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwasRunImputationMethodID, gwasRunKinshipID, gwasRunPopulationStructureID)
    insertedGwasRunID = insert_gwas_run(conn, args, gwas_run_obj)
    insertedGwasRunIDs.append(insertedGwasRunID)
  return insertedGwasRunIDs


def insert_gwas_result(conn, args, gwas_result):
  """Inserts gwas_result into database

  This function inserts a gwas_result into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    gwas_result (gwas_result): :ref:`gwas_result <gwas_result_class>` object
  
  Returns:
    int: GWAS result id
  """
  cur = conn.cursor()
  SQL = """INSERT INTO gwas_result(gwas_result_chromosome, basepair, gwas_result_gwas_run, pval, cofactor, _order, null_pval, model_added_pval, model, pcs)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_result_id;"""
  args = (gwas_result.c, gwas_result.b, gwas_result.r, gwas_result.p, gwas_result.o, gwas_result.d, gwas_result.n, gwas_result.a, gwas_result.m, gwas_result.s)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None


def insert_gwas_results_from_file(conn, args, speciesID, gwas_results_file, gwas_algorithm_ID, missing_snp_cutoff_value, missing_line_cutoff_value, imputationMethodID, genotypeVersionID, kinshipID, populationStructureID, minor_allele_frequency_cutoff_value):
  """Inserts a collection of GWAS results from a file into database

  This function inserts a collection of GWAS results from a file into a database

  Args:
    conn (psycopg2.extensions.connection): psycopg2 connection
    args (ArgumentParser namespace): user-defined arguments
    speciesID (int): :ref:`species_id <species_class>`
    gwas_results_file (str): absolute path to input file
    gwas_algorithm_ID (int): :ref:`gwas_algorithm_id <gwas_algorithm_class>`
    missing_snp_cutoff_value (numeric):
    missing_line_cutoff_value (numeric):
    imputationMethodID (int): :ref:`imputation_method_id <imputation_method_class>`
    genotypeVersionID (int): :ref:`genotype_version_id <genotype_version_class>`
    kinshipID (int): :ref:`kinship_id <kinship_class>`
    populationStructureID (int): :ref:`population_structure_id <population_structure_class>`
    minor_allele_frequency_cutoff_value (numeric):

  Returns:
    list of int: GWAS Result IDs
  """
  new_gwas_result_IDs = []
  df = pd.read_csv(gwas_results_file)
  for index, row in tqdm(df.iterrows(), total=df.shape[0], desc="GWAS Results"):
    trait = row['trait']
    traitID = find.find_trait(conn, args, trait)
    gwas_run_ID = find.find_gwas_run(conn, args, gwas_algorithm_ID, missing_snp_cutoff_value, missing_line_cutoff_value, imputationMethodID, traitID, row['nSNPs'], row['nLines'], genotypeVersionID, kinshipID, populationStructureID, minor_allele_frequency_cutoff_value)
    snp = row['SNP']
    snp_list = snp.split("_")
    chromosome = snp_list[0]
    chromosome = "chr"+str(chromosome)
    chromosomeID = find.find_chromosome(conn, args, chromosome, speciesID)
    basepair = snp_list[1]
    
    pcs = row['PCs']
    if type(pcs) == str:
      pcs_list = pcs.split(":")
      pcs_list = [int(x) for x in pcs_list]
    elif np.isnan(pcs):
      pcs_list = None
 
    new_gwas_result = gwas_result(chromosomeID, basepair, gwas_run_ID, row['pval'], row['cofactor'], row['order'], row['nullPval'], row['modelAddedPval'], row['model'], pcs_list)
    # if new_gwas_result:
    #   print("GWAS result entry successfully imported.")
    # else:
    #   print("GWAS result entry import FAILED.")
    new_gwas_result_ID = insert_gwas_result(conn, args, new_gwas_result)
    new_gwas_result_IDs.append(new_gwas_result_ID)
  return new_gwas_result_IDs
