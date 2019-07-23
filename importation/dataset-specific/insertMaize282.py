#!/usr/bin/env python
"""Hard-coded importation script for Maize282 Diverity Panel dataset.
The goal of this is to just get the data into the database and have a better grasp on how
to generalize importation given a strict file structure. This is a refactor of Molly Wohl's
original insertMaize282.py script.

The MLMM GWAS Pipeline (by Greg Ziegler) is the source for the majority of information.

GitHub: https://github.com/gziegler/runBatchMLMM.BellwetherGBMaize

There are two listed repos for this because it appears that the first, primary one was
removed.

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

Curated Values & Files:
  The following values must be known. The listed files are the best guesses. Any value marked with ⚠ is not yet confirmed. 
    * Species shortname: maize
    * Species binomial name: Zea mays
    * Species subspecies: None
    * Species variety: None
    * Population name: Maize282
    * Number of chromosomes: 10
    * Lines filename (.indv): <chromosome>_282_agpv4.012
    * Genotype version assembly name: B73 RefGen_v4
    * Genotype version annotation name: AGPv4
    * Reference genome line name (line name): 282set_B73
    * Phenotype filename(s): 5.mergedWeightNorm.LM.rankAvg.longFormat.csv
    * GWAS algorithm name: MLMM
    * Imputation method name: impute to major allele
    * Kinship algortihm name: van raden
    * Population structure algorithm name: Eigenstrat
    * Genotype filename(s) (.012): <chromosome>_282_agpv4.012
    * Variants filenames (.012.pos): <chromosome>_282_agpv4.012.pos
    * Kinship filename (.csv): 4.AstleBalding.synbreed.kinship.csv
    * Population structure filename (.csv): 4.Eigenstrat.population.structure.10PCs.csv
    * GWAS run filename (.csv): 9.mlmmResults.csv
    * GWAS results filename (.csv): 9.mlmmResults.csv
    * Missing SNP cutoff value: 0.2
    * Missing line cutoff value: 0.2
    * Minor allele frequency cutoff value: 0.1

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
  Lines:
  Genotype Version: < 1 sec
  Phenotypes: 
  Genotypes: 
  Variants: 
  Results:

  Estimated total time: 

"""

import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
from importation.util import find, insert
from importation.util.dbconnect import config, connect
from importation.util.models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result


if __name__ == '__main__':
  conn = connect()

  # ADD HARD-CODED VALUES FOR INDEPENDENT TABLES/OBJECTS

  # ADD LOCATIONS
  locations = []
  locations.append(location("United States", "Indiana", "West Lafayette", "PU"))
  locations.append(location("United States", "New York", None, "NY"))
  locations.append(location("United States", "Florida", None, "FL"))
  locations.append(location("United States", "Puerto Rico", None, "PR"))
  locations.append(location("United States", "North Carolina", None, "NC"))
  locations.append(location("South Africa", None, None, "SA"))
  locations.append(location("United States", "Missouri", None, "MO"))
  for place in locations:
    insert.insert_location(conn, place)
  # LOOK UP ID OF A HARD-CODED LOCATION USING find_chromosome()
  PUlocID = find.find_location(conn, 'PU')
  NYlocID = find.find_location(conn, "NY")
  FLlocID = find.find_location(conn, "FL")
  PRlocID = find.find_location(conn, "PR")
  NClocID = find.find_location(conn, "NC")
  SAlocID = find.find_location(conn, "SA")
  MOlocID = find.find_location(conn, "MO")

  # ADD A HARD-CODED SPECIES TO DB USING insert_species()
  soybeanSpecies = species('soybean', 'Glycine max', None, None)
  insertedSpeciesID = insert.insert_species(conn, soybeanSpecies)
  print("[ INSERT ]\t(%s)\t%s" % (insertedSpeciesID, str(soybeanSpecies)))
  mySpecies = species('maize', 'Zea mays', None, None)
  insertedSpeciesID = insert.insert_species(conn, mySpecies)
  print("[ INSERT ]\t(%s)\t%s" % (insertedSpeciesID, str(mySpecies)))
  maizeSpeciesID = find.find_species(conn, 'maize')
  print("[ FIND ]\t(%s)\t%s" % (maizeSpeciesID, '< species: maize >'))

  # ADD A HARD-CODED POPULATION TO DB USING insert_population()
  myPopulation = population('Maize282', maizeSpeciesID)
  insertedPopulationID = insert.insert_population(conn, myPopulation)
  print("[ INSERT ]\t(%s)\t%s" % (insertedPopulationID, str(myPopulation)))
  maize282popID = find.find_population(conn, 'Maize282')
  print("[ FIND ]\t(%s)\t%s" % (maize282popID, '< population: Maize282 >'))

  # ADD A HARD-CODED LINE TO DB USING insert_line()
  myLine = line(line_name='282set_B73', line_population=maize282popID)
  insertedLineID = insert.insert_line(conn, myLine)
  print("[ INSERT ]\t(%s)\t%s" % (insertedLineID, str(myLine)))
  B73lineID = find.find_line(conn, '282set_B73', maize282popID)
  print("[ FIND ]\t(%s)\t%s" % (B73lineID, '< line: Maize282 >'))

  # ADD NEW HARD-CODED GENOTYPE_VERSION TO DB
  myGenotypeVersion = genotype_version(genotype_version_name='B73 RefGen_v4_AGPv4_Maize282',
                                       genotype_version=315, reference_genome=B73lineID, genotype_version_population=maize282popID)
  B73_agpv4_maize282_versionID = insert.insert_genotype_version(conn, myGenotypeVersion)
  print("[ INSERT ]\t(%s)\t%s" % (B73_agpv4_maize282_versionID, str(myGenotypeVersion)))

  # ADD ALL CHROMOSOMES FOR A SPECIES TO DB
  insertedChromosomeIDs = insert.insert_all_chromosomes_for_species(conn, 10, maizeSpeciesID)
  print("[ INSERT ]\t%s\t%s" % (insertedChromosomeIDs, '\t10 (sID: %s)' % maizeSpeciesID))
  
  # GET LINES FROM SPECIFIED 012.indv FILE AND ADD TO DB
  insertedLineIDs = insert.insert_lines_from_file(conn, '../data/chr10_282_agpv4.012.indv', maize282popID)
  print("[ INSERT ]\t%s\t%s\t(pID:  %s)" % (insertedLineIDs, '../data/chr10_282_agpv4.012.indv', maize282popID))

  # GET VARIANTS FROM .012.pos FILE AND ADD TO  DB
  # Found the issue, the 'true' database on adriatic houses variants for ALL chromosomes
  # So, to fix that, we gotta loop through each chromosome file and add them
  # NOTE(timp): For when this is generalized to more than just Zea mays, there need to be a 
  # variable for the range instead because the number of chromosomes may differ between species
  for c in range(1, 11):
    chrShortname = 'chr' + str(c)
    chrId = find.find_chromosome(conn, chrShortname, maizeSpeciesID)
    filename = '../data/%s_282_agpv4.012.pos' % chrShortname
    # print("[ FIND ]\t(%s)\t%s" % (chrId, '< chromsome: %s >' % filename))
    insertedVariantIDs = insert.insert_variants_from_file(conn, filename, maizeSpeciesID, chrId)
    # print("num inserted variants:")
    # print(len(insertedVariantIDs))

  # ADD ALL GENOTYPES FROM A ONE-CHROMOSOME .012 FILE TO DB
  for c in range(1, 11):
    chrShortname = 'chr' + str(c)
    chrId = find.find_chromosome(conn, chrShortname, maizeSpeciesID)
    genoFilename = '../data/%s_282_agpv4.012' % chrShortname
    linesFilename = '../data/%s_282_agpv4.012.indv' % chrShortname
    # Example input file: chr1_282_agpv4.012.indv
    # 282set_33-16
    # 282set_38-11Goodman-Buckler
    # 282set_4226
    # 282set_4722
    # 282set_A188
    # 282set_A214NGoodman-Buckler
    # 282set_A239
    # 282set_A441-5
    # 282set_A554
    # ...
    # This is a list of all the lines that have been genotyped
    # AFAIK, this is 1:1 for the rows of each file, so row 1 of .indv contains the line of row 1 in .012

    insertedGenotypeIDs = insert.insert_genotypes_from_file(conn, genoFilename, linesFilename, chrId, maize282popID, B73lineID)
    # print("Inserted genotype IDs:")
    # print(insertedGenotypeIDs)
    # print("[ INSERT ]\t%s\t%s\t%s\t(cID: %s, pID: %s, lID: %s)" % (insertedGenotypeIDs, genoFilename, linesFilename, str(chrId), str(maize282popID), str(B73lineID)))

  # PARSE TRAITS FROM PHENOTYPE FILE AND ADD TO DB
  phenotypeRawData = pd.read_csv('../data/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', index_col=0)
  traits = list(phenotypeRawData)
  insertedTraitIDs = insert.insert_traits_from_traitlist(conn, traits)
  # print("num inserted traits:")
  # print(len(insertedTraitIDs))
  # print("Inserted trait IDs:")
  # print(insertedTraitIDs)
  
  # PARSE PHENOTYPES FROM FILE AND ADD TO DB
  # Example input file: 5.mergedWeightNorm.LM.rankAvg.longFormat.csv
  # Pedigree                      weight_FL06   weight_MO06   weight_NC06 ...
  # 282set_33-16                  299.8285      NA            247.08025
  # 282set_38-11Goodman-Buckler	  NA            157.62175     183.5531625
  # 282set_4226                   NA            NA            266.214
  # 282set_4722                   155.593625    130.501625    98.497
  # 282set_A188                   252.62675     255.4635      213.556125
  # 282set_A214NGoodman-Buckler	  NA            NA            202.21075
  # 282set_A239                   NA            225.50125     217.842
  # ...
  # It is a line for line listing of all the traits by year
  # This WILL be changed out for using phenotype (.ph) files instead

  insertedPhenoIDs = insert.insert_phenotypes_from_file(conn, '../data/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', maize282popID)
  # print("num phenotypes inserted:")
  # print(len(insertedPhenoIDs))
  # print("phenoIDs:")
  # print(insertedPhenoIDs)

  # ADD NEW HARD-CODED GROWOUT_TYPE TO DB
  greenhouse_GrowoutType = growout_type("greenhouse")
  greenhouse_GrowoutTypeID = insert.insert_growout_type(conn, greenhouse_GrowoutType)

  phenotyper_GrowoutType = growout_type("phenotyper")
  phenotyper_GrowoutTypeID = insert.insert_growout_type(conn, phenotyper_GrowoutType)

  field_GrowoutType = growout_type("field")
  field_GrowoutTypeID = insert.insert_growout_type(conn, field_GrowoutType)

  # LOOK UP ID OF A HARD-CODED GROWOUT_TYPE USING find_chromosome()
  fieldGrowoutTypeID = find.find_growout_type(conn, 'field')
  print("[ FIND ]\t(%s)\t%s" % (fieldGrowoutTypeID, '< growout_type: field >'))

  # ADD NEW HARD-CODED GROWOUT TO DB
  growouts = []
  growouts.append(growout("PU09", maize282popID, PUlocID, 2009, fieldGrowoutTypeID))
  growouts.append(growout("NY06", maize282popID, NYlocID, 2006, fieldGrowoutTypeID))
  growouts.append(growout("NY10", maize282popID, NYlocID, 2010, fieldGrowoutTypeID))
  growouts.append(growout("FL06", maize282popID, FLlocID, 2006, fieldGrowoutTypeID))
  growouts.append(growout("PR06", maize282popID, PRlocID, 2006, fieldGrowoutTypeID))
  growouts.append(growout("NC06", maize282popID, NClocID, 2006, fieldGrowoutTypeID))
  growouts.append(growout("PU10", maize282popID, PUlocID, 2010, fieldGrowoutTypeID))
  growouts.append(growout("SA06", maize282popID, SAlocID, 2006, fieldGrowoutTypeID))
  growouts.append(growout("MO06", maize282popID, MOlocID, 2006, fieldGrowoutTypeID))
  insertedGrowoutIDs = []
  for growout in growouts:
    print("-------------\t%s" % str(growout))
    insertedGrowoutIDs.append(insert.insert_growout(conn, growout))
  print("[ INSERT ]\t%s\t(new growout)" % (insertedGenotypeIDs) )
  
  # ADD NEW HARD-CODED GWAS_ALGORITHM TO DB
  gwasAlgorithms = []
  gwasAlgorithms.append(gwas_algorithm("MLMM"))
  gwasAlgorithms.append(gwas_algorithm("EMMAx"))
  gwasAlgorithms.append(gwas_algorithm("GAPIT"))
  gwasAlgorithms.append(gwas_algorithm("FarmCPU"))
  newGWASalgorithmIDs = []
  for algorithm in gwasAlgorithms:
    newGWASalgorithmIDs.append(insert.insert_gwas_algorithm(conn, algorithm))
  print("[ INSERT ]\t%s\t(new gwas algorithm IDs)" % (newGWASalgorithmIDs) )
  newGWASalgorithm = find.find_gwas_algorithm(conn, 'MLMM')


  # ADD NEW HARD-CODED IMPUTATION_METHOD TO DB
  newImputationMethods = []
  newImputationMethods.append(imputation_method("impute to major allele"))
  newImputationMethods.append(imputation_method("impute to minor allele"))
  newImputationMethods.append(imputation_method("impute to average allele"))
  newImputationMethods.append(imputation_method("IMPUTE"))
  newImputationMethods.append(imputation_method("BEAGLE"))
  for im in newImputationMethods:
    insert.insert_imputation_method(conn, im)
  
  # ADD NEW HARD-CODED KINSHIP_ALGORITHM TO DB
  kinshipAlgorithms = []
  kinshipAlgorithms.append(kinship_algorithm("loiselle"))
  kinshipAlgorithms.append(kinship_algorithm("van raden"))
  kinshipAlgorithms.append(kinship_algorithm("Synbreed_realizedAB"))
  newKinshipAlgorithmIDs = []
  for algorithm in kinshipAlgorithms:
    newKinshipAlgorithmIDs.append(
        insert.insert_kinship_algorithm(conn, algorithm))
  print("[ INSERT ]\t%s\t(new kinship algorithm IDs)" % (newKinshipAlgorithmIDs))
  # LOOK UP ID OF A HARD-CODED KINSHIP_ALGORITHM USING find_kinship_algorithm()
  VanRadenID = find.find_kinship_algorithm(conn, "van raden")
  print("Van Raden kinship alg ID:")
  print(VanRadenID)  

  # ADD NEW HARD-CODED KINSHIP TO DB
  newKinship = kinship(VanRadenID, "../data/4.AstleBalding.synbreed.kinship.csv")
  newKinshipID = insert.insert_kinship(conn, newKinship)
  print("New kinship ID:")
  print(newKinshipID)

  # ADD NEW HARD-CODED POPULATION_STRUCTURE_ALGORITHM TO DB
  newPopulationStructures = []
  newPopulationStructures.append(population_structure_algorithm("Eigenstrat"))
  newPopulationStructures.append(population_structure_algorithm("STRUCTURE"))
  newPopulationStructures.append(population_structure_algorithm("FastSTRUCTURE"))
  for ps in newPopulationStructures:
    insert.insert_population_structure_algorithm(conn, ps)

  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE_ALGORITHM USING find_population_structure_algorithm()
  EigenstratID = find.find_population_structure_algorithm(conn, "Eigenstrat")
  print("Eigenstrat pop str alg ID:")
  print(EigenstratID)

  # ADD NEW HARD-CODED POPULATION_STRUCTURE TO DB
  # Example input file: 4.Eingenstrat.population.structure.10PCs.csv
    # Line                         	      V1          	 V2	          V3 ...
    # 282set_4226                   -0.002298602  -0.029693879   0.008527265
    # 282set_4722                   -0.003785163	-0.083527265	-0.059586105
    # 282set_33-16                   0.000222197	-0.035755785   0.017007817
    # 282set_38-11Goodman-Buckler   -0.026698262	-0.053115302	-0.01159794
    # 282set_A188                    0.002520617	-0.041387288	-0.011656126
    # 282set_A239                   -0.024217977	-0.038008255   0.033222018
    # ...
    # The number of columns is one more than the number of PCs in filename
  newPopulationStructure = population_structure(EigenstratID, "../data/4.Eigenstrat.population.structure.10PCs.csv")
  newPopulationStructureID = insert.insert_population_structure(conn, newPopulationStructure)
  print("New population structure ID:")
  print(newPopulationStructureID)

  # LOOK UP ID OF A HARD-CODED GWAS_ALGORITHM
  MLMMalgorithmID = find.find_gwas_algorithm(conn, "MLMM")
  print("MLMM algorithm ID:")
  print(MLMMalgorithmID)

  # LOOK UP ID OF A HARD-CODED GENOTYPE_VERSION
  B73_agpv4_maize282_versionID = find.find_genotype_version(conn, "B73 RefGen_v4_AGPv4_Maize282")
  print("B73 agpv4 maize282 genotype version: ")
  print(B73_agpv4_maize282_versionID)  

  # LOOK UP ID OF A HARD-CODED IMPUTATION_METHOD
  majorAlleleImputationID = find.find_imputation_method(conn, "impute to major allele")
  print("major allele imputation ID: ")
  print(majorAlleleImputationID)  

  # LOOK UP ID OF A HARD-CODED KINSHIP
  # NOTE(timp): I could not find this file, but I found a R data file (.rda) that may contain the information.
  #             Although, the data may not be in the correct format.
  #             The temporary file is the one with 'export' in its name.
  # kinshipID = find.find_kinship(conn, "/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv")
  kinshipID = find.find_kinship(conn, "../data/4.AstleBalding.synbreed.kinship.csv")
  print("kinshipID: ")
  print(kinshipID)  

  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE
  populationStructureID = find.find_population_structure(conn, "../data/4.Eigenstrat.population.structure.10PCs.csv")
  print("population structure ID: ")
  print(populationStructureID)

  # PARSE GWAS_RUNS FROM FILE AND ADD TO DB
  # NOTE(timp): Could not find file or possible equivalent
  insertedGwasRunIDs = insert.insert_gwas_runs_from_gwas_results_file(conn,
                                                                      '../data/9.mlmmResults.csv',
                                                                      MLMMalgorithmID,
                                                                      B73_agpv4_maize282_versionID,
                                                                      0.2,
                                                                      0.2,
                                                                      0.1,
                                                                      majorAlleleImputationID,
                                                                      kinshipID,
                                                                      populationStructureID)
  print("Inserted gwas_run IDs:")
  print(insertedGwasRunIDs)

  # PARSE GWAS_RESULTS FROM FILE AND ADD TO DB
  # NOTE(timp): Could not find file or possible equivalent
  insertedGwasResultIDs = insert.insert_gwas_results_from_file(conn,
                                                              maizeSpeciesID,
                                                              '../data/9.mlmmResults.csv',
                                                              MLMMalgorithmID,
                                                              0.2,
                                                              0.2,
                                                              majorAlleleImputationID,
                                                              B73_agpv4_maize282_versionID,
                                                              kinshipID,
                                                              populationStructureID,
                                                              0.1)
  print("Inserted gwas result IDs: ")
  print(insertedGwasResultIDs)
