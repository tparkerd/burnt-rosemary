#########################
Maize 282 Diversity Panel
#########################

*uploaded: 9999-12-31*

Description
  *todo*


***********
Known Files
***********

.012
=========
Contains a list of all the genotypes (allele calls) for each line given a
chromosome.

Usage: ``insert.insert_genotypes_from_file()``

- chr1_282_agpv4.012
- chr2_282_agpv4.012
- chr3_282_agpv4.012
- chr4_282_agpv4.012
- chr5_282_agpv4.012
- chr6_282_agpv4.012
- chr7_282_agpv4.012
- chr8_282_agpv4.012
- chr9_282_agpv4.012
- chr10_282_agpv4.012

Example file:

.. code-block:: bash
  :linenos:

  # chr10_282_agpv4.012 (277 rows, 422,653 columns)
  0 2 2 0 2 2 2 2 0 0 2 0 ...
  1 2 2 0 2 2 2 2 2 0 2 0 ...
  2 2 2 0 2 2 2 2 0 0 2 0 ...
  3 0 0 0 0 0 0 0 1 0 1 0 ...
  4 1 1 0 1 2 2 2 0 0 2 0 ...
  5 0 0 0 0 0 0 0 0 0 0 0 ...
  ...

.012.indv
=========

Contains a list of all the lines given a chromosome

Usage: ``insert.insert_genotypes_from_file()``

- chr1_282_agpv4.012.indv
- chr2_282_agpv4.012.indv
- chr3_282_agpv4.012.indv
- chr4_282_agpv4.012.indv
- chr5_282_agpv4.012.indv
- chr6_282_agpv4.012.indv
- chr7_282_agpv4.012.indv
- chr8_282_agpv4.012.indv
- chr9_282_agpv4.012.indv
- chr10_282_agpv4.012.indv

Example file:

.. code-block:: bash
  :linenos:

  # chr10_282_agpv4.012.indv (277 rows)
  282set_33-16
  282set_38-11Goodman-Buckler
  282set_4226
  282set_4722
  282set_A188
  282set_A214NGoodman-Buckler
  282set_A239
  282set_A441-5
  282set_A554
  282set_A556
  282set_A6
  282set_A619
  282set_A632
  ...

.012.pos
=========

Lists all the SNP positions for each chromosome. Each row of this file maps to
a column of its corresponding ``.012`` file.

Usage: ``insert.insert_variants_from_file()``

- chr1_282_agpv4.012.pos
- chr2_282_agpv4.012.pos
- chr3_282_agpv4.012.pos
- chr4_282_agpv4.012.pos
- chr5_282_agpv4.012.pos
- chr6_282_agpv4.012.pos
- chr7_282_agpv4.012.pos
- chr8_282_agpv4.012.pos
- chr9_282_agpv4.012.pos
- chr10_282_agpv4.012.pos

Example File:

.. code-block:: bash
  :linenos:

  # chr10_282_agpv4.012.pos (422,652 row)
  10	94641
  10	94642
  10	94653
  10	94668
  10	94901
  10	104377
  10	104384
  10	108473
  10	108480
  10	108494


5.mergedWeightNorm.LM.rankAvg.longFormat.csv
--------------------------------------------

Phenotype data
