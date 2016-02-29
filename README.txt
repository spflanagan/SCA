SCA Analysis README

This file contains information about each program


*******************************************process_alleles_files**************************************************
The purpose of process_alleles_files is to process the alleles files output from stacks (**.alleles.tsv) and reformat them to be tab-delimited with the columns IndID, LocusID, Haplotype1, Count1, Haplotype2. If an individual has more than two allele calls for that locus the Haplotype and Count will be added on as additional columns in that locus' row.

This file can be run interactively or can be run at the command line with flags. You must give it the input alleles filename with the path and the output filename with the path. An optional input is a whitelist file with a list of Catalog IDs (one ID per row).

**********************************************convert_matches*******************************************************
The purpose of convert_matches is to process the matches files output from stacks (**.matches.tsv) and reformat them to be tab-delimited files with the following columns: CatalogID, Allele1, Allele1 Count, Allele2, Allele2 Count.

This file can be run interactively or can be run at the command line with flags. You must give it the input alleles filename with the path and the output filename with the path. An optional input is a whitelist file with a list of Catalog IDs (one ID per row).