#######################
Setaria Diversity Panel
#######################

The Setaria Diversity Panel


************
Known Issues
************

During the analysis process, two lines could not be mapped from the original versions that used a *XXXX_XXXX_XXXX_XXXX* format to more human friendly named such as *TB_0003*. The two specific lines were ``ACPAW_ACPAW_ACPAW_ACPAW`` and ``ACPAX_ACPAX_ACPAX_ACPAX``. The human friendly names are ``TB_0048`` or ``TB_0261``. We cannot determine which one maps to which. Therefore, they were exluded from the analysis. The ``REF_REF_REF_REF`` line could not be found, but according to Greg Zeigler, the reference line is ``B100``.


Genotype, Individual, and Position (.012, .indv, .pos) files were generated
from the original files:

``2.from12.setaria.maf0.1.maxMissing0.1.allLines.012``
``2.from12.setaria.maf0.1.maxMissing0.1.allLines.012.indv``
``2.from12.setaria.maf0.1.maxMissing0.1.allLines.012.pos``

The population structure file (6.Eigenstrat.population.structure.50PCs.csv)
has been modified from the original copy because the column names were
shifted left by one. Therefore, I added the 'Pedigree' column title and
shifted the V# headings to the right, adding V50 to the final column.

The mappings for weird namings for line/pedigrees is in the file
Setaria_597_diversity_samples.csv
Still need to include the conversion in the script. Or! Modify the existing
files to use them expected pedigree labels instead of the strange
XXXX_XXXX_XXXX_XXXX construction.



*****
Files
*****

List of Files
-------------

:Genotypes: chr{1..9}_setaria.012
:Lines: chr{1..9}_setaria.012.indv
:Phenotypes / Traits: 2.Setaria_IR_2016_datsetset_GWAS.BLUPsandBLUEs.csv
:Kinship: 6.AstleBalding.synbreed.kinship.csv'
:Population Structure: 6.Eigenstrat.population.structure.50PCs.csv
:GWAS Run / Results: 11.allStomataresults.csv

.. literalinclude:: ../../../../data/setaria.json
    :caption: :download:`Download setaria.json <../../../../data/setaria.json>`
    :name: setaria.json
    :language: json
    :linenos:

.012 (Genotype file)
Format: