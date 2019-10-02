#!/usr/bin/env python
""""""
import os
import logging
import argparse
from datetime import datetime as dt

class Experiment():
  def __init__(self):
    pass

def validate(args):
  try:
    # Open configuration file
    with open(args.filename, 'r') as ifp:
      cfg = json.load(ifp)
      required_fields = [ "species_shortname",
                          "species_binomial_name",
                          "species_subspecies",
                          "species_variety",
                          "population_name",
                          "number_of_chromosomes",
                          "genotype_version_assembly_name",
                          "genotype_version_annotation_name",
                          "reference_genome_line_name",
                          "phenotype_filename" ]
      
      missing_keys = []
      for k in required_fields:
        if k not in cfg:
          missing_keys.append(k)
      if missing_keys:
        raise KeyError(f'The following keys are required. Please include them in your JSON configuration: {missing_keys}')

      logging.info(f"Configuration file is valid. Verifying that all files exist.")

def design(args):
  # =======================================
  # ========== Experiment Design ==========
  # =======================================
  # What the database needs in order to create an 'experiment' is the following
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
  species_shortname = cfg['species_shortname']
  species_binomial = cfg['species_binomial_name']
  species_subspecies = cfg['species_subspecies']
  species_variety = cfg['species_variety']
  # Population
  population_name = cfg['population_name']
  # Chromosome
  chromosome_count = cfg['number_of_chromosomes']
  # Line
  # NOTE(tparker): Can use any chromosome, as they are the same for each.
  # In the future, this the extraneous copies of the lines may be removed
  # and there will be one specific line file, much like the phenotype files
  lines_filename = Template('${cwd}/${chr}_${shortname}.012.indv')
  # # Genotype Version
  # # NOTE(tparker): This is possibly just the info about the reference genome
  # #                It is likely included with the VCF genotype file (.012).
  # genotype_version_assembly_name = cfg['genotype_version_assembly_name']
  # genotype_version_annotation_name = cfg['genotype_version_annotation_name']
  # reference_genome_line_name = cfg['reference_genome_line_name']
  # # Growout, Type, and Location
  # # NOTE(tparker): Unknown at this time
  # ## Location
  # ## Type
  # ## Growout
  # #
  # # Traits
  # # Allow for more than on phenotype files
  # if isinstance(cfg["phenotype_filename"], list):
  #   phenotype_filenames = [ f'{args.wd}/{filename}' for filename in cfg['phenotype_filename'] ]
  # else:
  #   phenotype_filenames = [ f'{args.wd}/{configuration["phenotype_filename"]}']

  # # Model Construction & Insertion
  # # Species
  # s = species(species_shortname, species_binomial, species_subspecies, species_variety)
  # species_id = insert.insert_species(conn, args, s)
  # logging.debug(f'[Insert]\tSpecies ID\t{species_id}, {s}')
  # # Population
  # p = population(population_name, species_id)
  # population_id = insert.insert_population(conn, args, p)
  # logging.debug(f'[Insert]\tPopulation ID\t{population_id}: {p}')
  # # Chromosome
  # chromosome_ids = insert.insert_all_chromosomes_for_species(conn, args, chromosome_count, species_id)
  # logging.debug(f'[Insert]\tChromosome IDs\t{chromosome_ids}')
  # # Line
  # working_filepath = lines_filename.substitute(dict(chr="chr1", cwd=f"{args.wd}", shortname=species_shortname))
  # try:
  #   if not os.path.isfile(working_filepath):
  #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), working_filepath)
  # except:
  #   raise

  # line_ids = insert.insert_lines_from_file(conn, args, working_filepath, population_id) # hard-coded substitue until just one file is used for lines
  # logging.debug(f'[Insert]\tLine IDs\t{line_ids}')
  # # Genotype Version
  # reference_genome_id = find.find_line(conn, args, reference_genome_line_name, population_id)
  # logging.debug(f'[Insert]\tReference Genome ID\t{reference_genome_id}, ({reference_genome_line_name}, {population_id})')
  # gv = genotype_version(genotype_version_name = genotype_version_assembly_name,
  #                       genotype_version = genotype_version_annotation_name,
  #                       reference_genome = reference_genome_id,
  #                       genotype_version_population = population_id)
  # genotype_version_id = insert.insert_genotype_version(conn, args, gv)
  # logging.debug(f'[Insert]\tGenome Version ID\t{genotype_version_id}')
  # if genotype_version_id is None:
  #   raise Exception(f'Genotype version is None for parameters: {pformat(gv)}')
  
  # # Growout, Type, and Location
  # # NOTE(tparker): Unknown at this time
  # ## Location
  # ## Type
  # ## Growout

  # # Traits
  # # Go through all the phenotype files available for the dataset and insert
  # # the recorded traits for each.
  # for phenotype_filepath in phenotype_filenames:
  #   try:
  #     if not os.path.isfile(phenotype_filepath):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), phenotype_filepath)
  #   except:
  #     raise
  #   traits = list(pd.read_csv(phenotype_filepath, index_col=0))
  #   trait_ids = insert.insert_traits_from_traitlist(conn, args, traits, phenotype_filepath)
  #   logging.debug(f'[Insert]\tTrait IDs for {phenotype_filepath}\t{trait_ids}')

def collect(args):
  pass
  # # ===========================================
  # # ========== Experiment Collection ==========
  # # ===========================================
  # # Phenotype (external source?)
  # #       This needs to be standardized to a .pheno filetype.
  # #       For now, it is the longFormat for the Maize282 datset
  # #       5.mergedWeightNorm.LM.rankAvg.longFormat.csv, but for Setaria will be 
  # # Genotype (VCF output)
  # # Variant (VCF output)

  # # Expected User Input
  # # Phenotype
  # # NOTE(tparker): Define in earlier stage
  # # Genotype
  # genotype_filename = Template('${cwd}/${chr}_${shortname}.012')
  # # Variants
  # variants_filename = Template('${cwd}/${chr}_${shortname}.012.pos')

  # # Model Construction & Insertion
  # # Phenotype
  # for phenotype_filepath in phenotype_filenames:
  #   try:
  #     if not os.path.isfile(phenotype_filepath):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), phenotype_filepath)
  #   except:
  #     raise
  #   phenotype_ids = insert.insert_phenotypes_from_file(conn, args, phenotype_filepath, population_id, phenotype_filepath)
  #   logging.debug(f'[Insert]\tPhenotype IDs for {phenotype_filepath}\t{phenotype_ids}')

  # # Genotype
  # for c in range(1, chromosome_count + 1):
  #   chromosome_shortname = 'chr' + str(c)
  #   chromosome_id = find.find_chromosome(conn, args, chromosome_shortname, species_id)
  #   geno_filename = genotype_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.wd}', shortname=species_shortname))
  #   line_filename = lines_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.wd}', shortname=species_shortname))
  #   try:
  #     if not os.path.isfile(geno_filename):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), geno_filename)
  #     if not os.path.isfile(line_filename):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), line_filename)
  #   except:
  #     raise
  #   genotype_ids = insert.insert_genotypes_from_file(conn = conn,
  #                                                   args = args,
  #                                                   genotypeFile = geno_filename,
  #                                                   lineFile = line_filename,
  #                                                   chromosomeID = chromosome_id,
  #                                                   populationID = population_id,
  #                                                   genotype_versionID = genotype_version_id
  #                                                   )
  # # Variants
  # for c in range(1, chromosome_count + 1):
  #   chromosome_shortname = 'chr' + str(c)
  #   chromosome_id = find.find_chromosome(conn, args, chromosome_shortname, species_id)
  #   variant_filename = variants_filename.substitute(dict(chr=chromosome_shortname, cwd=f'{args.wd}', shortname=species_shortname))
  #   try:
  #     if not os.path.isfile(variant_filename):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), variant_filename)
  #   except:
  #     raise
  #   # insert.insert_variants_from_file(conn,
  #   #                                  args,
  #   #                                  variant_filename,
  #   #                                  species_id,
  #   #                                  chromosome_id)

  #   # NOTE(tparker): Changed variant insertion to the async version
  #   insert.insert_variants_from_file_async(conn,
  #                                         args,
  #                                         variant_filename,
  #                                         species_id,
  #                                         chromosome_id)

def parse_options():
  """Function to parse user-provided options from terminal"""
  parser = argparse.ArgumentParser(description="Template sub-module of gwas data importer")
  parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("-V", "--version", action="version", version='%(prog)s 1.0.0')
  parser.add_argument("-f", "--force", action="store_true", default=False, help="Force file creation. Overwrite any existing files.")
  parser.add_argument("--log", action="store", default=f"{dt.today().strftime('%Y-%m-%d')}-{os.path.basename(__file__)}.log", help="")
  # parser.add_argument('files', metavar='FILES', type=str, nargs='+', help='List of .nsihdr files')
  args = parser.parse_args()

  # Configure logging, stderr and file logs
  logging_level = logging.INFO
  if args.verbose:
    logging_level = logging.DEBUG

  logFormatter = logging.Formatter("%(asctime)s - [%(levelname)-4.8s] - %(filename)s %(lineno)d - %(message)s")
  rootLogger = logging.getLogger()
  rootLogger.setLevel(logging_level)

  fileHandler = logging.FileHandler(args.log)
  fileHandler.setFormatter(logFormatter)
  rootLogger.addHandler(fileHandler)

  consoleHandler = logging.StreamHandler()
  consoleHandler.setFormatter(logFormatter)
  rootLogger.addHandler(consoleHandler)
  return args

if __name__ == "__main__":
  args = parse_options()
  try:
    validate(args)
    design(args)
    collect(args)
  except Exception as err:
    logging.error(f'{type(err).__name__} - {err}'.replace('\n\n', '\t'))
    
  
  