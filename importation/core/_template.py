#!/usr/bin/env python
""""""
import os
import logging
import argparse
from datetime import datetime as dt

def validate(args):
  pass

def design(args):
  pass

def collect(args):
  pass

def parse_options():
  """Function to parse user-provided options from terminal"""
  parser = argparse.ArgumentParser(description="Template sub-module of gwas data importer")
  parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
  parser.add_argument("-V", "--version", action="version", version='%(prog)s 1.0.0')
  parser.add_argument("-f", "--force", action="store_true", default=False, help="Force file creation. Overwrite any existing files.")
  parser.add_argument("--log", action="store", default=f"{dt.today().strftime('%Y-%m-%d')}-{os.path.basename(__file__)}", help="")
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
    
  
  