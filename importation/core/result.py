#!/usr/bin/env python
""""""
import os
import logging
import argparse

def validate():
  pass

def design():
  pass

def collect():
  pass

  # # =============================================
  # # ================== Results ==================
  # # =============================================
  # # GWAS Run
  # # GWAS Results

  # # Expected User Input
  # # GWAS Run & results
  # if isinstance(dp['gwas_results_filename'], list):
  #   gwas_filenames = [ f'{args.wd}/{filename}' for filename in dp['gwas_results_filename'] ] # allows for more than one gwas results/run file
  # else:
  #   gwas_filenames = [ f'{args.wd}/{dp["gwas_results_filename"]}' ]
  # # The following values (0.2, 0.2, and 0.1) were all taken from the Maize282 import
  # # NOTE(tparker): Make sure to double check with Greg on what the true values should be
  # #                Also, double check the source of the pipeline to see if there is any
  # #                indication what the values shoudl be.
  # missing_snp_cutoff_value = dp['missing_SNP_cutoff_value']
  # missing_line_cutoff_value = dp['missing_line_cutoff_value']
  # minor_allele_frequency_cutoff_value = dp['minor_allele_frequency_cutoff_value']

  # # Model Construction & Insertion
  # # GWAS Run
  # # NOTE(tparker): Check with Greg on what the imputation method was used. I believe it was
  # #                set by someone named Sujan because imputation was done beforehand
  # for gwas_filename in gwas_filenames:  
  #   try:
  #     if not os.path.isfile(gwas_filename):
  #       raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), gwas_filename)
  #   except:
  #     raise
  #   imputation_method_id = find.find_imputation_method(conn, args, imputation_method_name)
  #   gwas_run_ids = insert.insert_gwas_runs_from_gwas_results_file(conn,
  #                                                                 args,
  #                                                                 gwas_filename,
  #                                                                 gwas_algorithm_id,
  #                                                                 genotype_version_id,
  #                                                                 missing_snp_cutoff_value,
  #                                                                 missing_line_cutoff_value,
  #                                                                 minor_allele_frequency_cutoff_value,
  #                                                                 imputation_method_id,
  #                                                                 kinship_id,
  #                                                                 population_structure_id)
  #   # GWAS Results
  #   gwas_result_ids = insert.insert_gwas_results_from_file(conn = conn,
  #                                                          args = args,
  #                                                          speciesID = species_id,
  #                                                          gwas_results_file = gwas_filename,
  #                                                          gwas_algorithm_ID = gwas_algorithm_id,
  #                                                          missing_snp_cutoff_value = missing_snp_cutoff_value,
  #                                                          missing_line_cutoff_value = missing_line_cutoff_value,
  #                                                          imputationMethodID = imputation_method_id,
  #                                                          genotypeVersionID = genotype_version_id,
  #                                                          kinshipID = kinship_id,
  #                                                          populationStructureID = population_structure_id,
  #                                                          minor_allele_frequency_cutoff_value = minor_allele_frequency_cutoff_value)


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
  
  