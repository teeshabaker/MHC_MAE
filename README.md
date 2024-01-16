## SUMMARY
This repository contains the R code used for my PhD thesis, "Improving Detection and Quantification of Major Histocompatibility Complex (MHC)-Presented Immunopeptides for Vaccine Development". 

Following code files are available in this repository.

1. MAE functions.R - This file contains all the helper functions that are used throughout the codebase
2. MAE_0B_fragpipe_PCA_exploration.rmd - This file is for running a principal component analysis (PCA) for data exploration
3. MAE_1A_byonic_processing.rmd - This file is for processing the excel output files from Byonic
4. MAE_1B_fragpipe_processing.rmd - This file is for processing the excel output files from Fragpipe
5. MAE_2_netMHCpan_processing.rmd - This file takes netMHCpan files as input, and creates a summarized list of binders, and their corresponding allele, based on binding strength as output
6. MAE_3B_epitope quantification.rmd - This file is for creating a quantification table for all salmonella epitopes identified across the multiple exposure conditions
7. MAE_4_GibbsCluster EPS to PDF to PNG.rmd - This file is for creating PNG files from the EPS files created by GibbCluster
8. MAE_5_combining GibbsCluster core files to create common motifs.rmd - This file is for combining sequences from mulitple .core files from GibbsCluster to create a summarized FASTA file
9. MAE_6_combining FIMO IDd similar sequences from protein of interest.rmd - This file is used to combine FIMO outputs for all common binding motifs for a single protein of interest and creates a FASTA file and peptie list file for visualizing sequence coverage with MStools
10. MAE_7_MATRIX2MEME.rmd - This file is for using a combined .mat frequency matrix from MAE_3 as input to create a .txt frequency matrix as the output that is readable by MATRIX2MEME from MEMESuite
11. items_to_remove_for_1000_genome.rmd - This file contains the item list that are removed prior to the analysis

## PREREQUISITES

HARDWARE:

- Good CPU (i5+), multicore preferred if we want to use parallel processing for speed
- Minimum 16GB RAM
- Available storage (10GB+)

SOFTWARE:

- R (version 4.0+)
- R Studio (version 2023+)
- Various R libraries as specified in individual code files


## CITING

If you use this code in your research, please use the following citation for my PhD thesis:

Baker, Teesha "Improving Detection and Quantification of Major Histocompatibility Complex (MHC)-Presented Immunopeptides for Vaccine Development", PhD thesis, University of British Columbia, February. 2024.


## LICENSE

The code is distributed under the GPL license (v3), please see the license file for more details.

## DISCLAIMER: USE OF CODE FROM THIS REPOSITORY

The code contained in this repository is provided for informational and educational purposes only. Users are hereby informed that the code has been sanitized and cleaned for public dissemination, and it is not an exact representation of the code utilized in the experiments conducted for the associated thesis.

The purpose of sharing this code is to offer insights, examples, and reference material. Users should exercise caution and due diligence when employing or relying on any code from this repository. The repository's content is not intended for direct application in production environments or critical systems without careful review, modification, and testing.

The author and contributors of this repository make no representations or warranties regarding the accuracy, completeness, or fitness for a particular purpose of the code provided. Any use of the code is at the user's own risk, and the author and contributors shall not be liable for any damages, losses, or adverse consequences arising out of or in connection with the use of the code.

Furthermore, users are encouraged to refer to the original thesis and associated documentation for a comprehensive understanding of the experiments, methodologies, and findings. The cleaned code in this repository is not a substitute for the detailed information provided in the academic work.

By accessing and using the code from this repository, users acknowledge and agree to the terms of this disclaimer
