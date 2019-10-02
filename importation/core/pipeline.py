#!/usr/bin/env python
""""""
import os
import logging
import argparse

def validate():
  pass

def design():
  pass

 
  # # # =====================================
  # # # ========== Pipeline Design ==========
  # # # =====================================
  # # # GWAS Algorithm: "MLMM", "EMMAx", "GAPIT", "FarmCPU"
  # # # Imputation Method: "impute to major allele"
  # # # Kinship Algorithm: "loiselle"
  # # # Population Structure Algorithm: "Eigenstrat"

  # # Expected User Input
  # # GWAS Algorithm
  # gwas_algorithm_name = dp['gwas_algorithm_name'] # According to Greg's README
  # # Imputation Method
  # imputation_method_name = dp['imputation_method_name'] # Unknown, apparently it was done by someone named Sujan
  # # Kinship Algorithm
  # kinship_algorithm_name = dp['kinship_algortihm_name'] # Placeholder, I don't know the exact string that should be used
  # # Population Structure Algorithm
  # population_structure_algorithm_name = dp['population_structure_algorithm_name'] # This is a guess based on filename

  # # Model Construction & Insertion
  # # GWAS Algorithm
  # ga = gwas_algorithm(gwas_algorithm_name)
  # gwas_algorithm_id = insert.insert_gwas_algorithm(conn, args, ga)
  # # Imputation Method
  # im = imputation_method(imputation_method_name)
  # imputation_method_id = insert.insert_imputation_method(conn, args, im)
  # # Kinship Algorithm
  # ka = kinship_algorithm(kinship_algorithm_name)
  # kinship_algorithm_id = insert.insert_kinship_algorithm(conn, args, ka)
  # # Population Structure Algorithm
  # psa = population_structure_algorithm(population_structure_algorithm_name)
  # population_structure_algorithm_id = insert.insert_population_structure_algorithm(conn, args, psa)


def collect():
  pass

  # # =========================================
  # # ========== Pipeline Collection ==========
  # # =========================================
  # # Kinship
  # # Setaria Kinship is stored in:
  # ## /shares/ibaxter_share/gziegler/SetariaGWASPipeline/data/genotype/6.AstleBalding.synbreed.kinship.rda 
  # ## Exported the file to CSV using R
  # ### load('6.AstleBalding.synbreed.kinship.rda')
  # ### write.csv(kinship, '6.AstleBalding.synbreed.kinship.csv')
  # # Population Structure

  # # Expected User Input
  # # Kinship
  # # NOTE(tparker): Currently the database just stores the filename.
  # #                There is no service to upload the file to database's
  # #                host, so there's no single location to find the file
  # #                I would like to find out why this is the case and if 
  # #                it would just be better to store it in the database and
  # #                allow the user to export the table themselves as a CSV.
  # kinship_filepath = f'{args.wd}/{dp["kinship_filename"]}'
  # # Population Structure
  # # NOTE(tparker): Same reasoning as the kinship file. There should be a way 
  # #                for the data to be stored in the database, not a 
  # population_structure_filepath = f'{args.wd}/{dp["population_structure_filename"]}'

  # # Model Construction & Insertion
  # # Kinship
  # try:
  #   if not os.path.isfile(kinship_filepath):
  #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), kinship_filepath)
  #   if not os.path.isfile(population_structure_filepath):
  #     raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), population_structure_filepath)
  # except:
  #   raise
  # k = kinship(kinship_algorithm_id, kinship_filepath)
  # kinship_id = insert.insert_kinship(conn, args, k)
  # # Population Structure
  # ps = population_structure(population_structure_algorithm_id, population_structure_filepath)
  # population_structure_id = insert.insert_population_structure(conn, args, ps)



def parse_options():
  """Function to parse user-provided options from terminal"""
  parser = argparse.ArgumentParser(description="This tool converts a NSI project from 32-bit float to 16-bit unsigned integer format, and it extracts the midslice and generates a side-view projection of the volume.")
  parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("-V", "--version", action="version", version='%(prog)s 1.0.0')
  parser.add_argument("-f", "--force", action="store_true", default=False, help="Force file creation. Overwrite any existing files.")
  parser.add_argument("--log", action="store", default=f"{dt.today().strftime("%Y-%m-%d")-{os.path.basename(__file__)}}")
  # parser.add_argument('files', metavar='FILES', type=str, nargs='+', help='List of .nsihdr files')
  args = parser.parse_args()

  # Configure logging, stderr and file logs
  logging_level = logging.INFO
  if args.verbose:
    logging_level = logging.DEBUG

  logFormatter = logging.Formatter("%(asctime)s - [%(levelname)-4.8s]  %(message)s")
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
  
  