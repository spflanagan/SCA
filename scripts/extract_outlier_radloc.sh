#!/bin/bash

#run from scripts/ directory
#~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
#	-c ../results/stacks/batch_1.catalog.tags.tsv \
#	-l ../results/biallelic_outliers/AJ_outliers.txt \
#	-o ../results/biallelic_outliers/rad_locus/aj_outlier_loci.fasta 
	
~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
	-c ../results/stacks/batch_1.catalog.tags.tsv \
	-l ../results/biallelic_outliers/FM_1outliers.txt \
	-o ../results/biallelic_outliers/rad_locus/fm_1outlier_loci.fasta 

~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
	-c ../results/stacks/batch_1.catalog.tags.tsv \
	-l ../results/biallelic_outliers/MO_1outliers.txt \
	-o ../results/biallelic_outliers/rad_locus/mo_1outlier_loci.fasta 
	
~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
	-c ../results/stacks/batch_1.catalog.tags.tsv \
	-l ../results/biallelic_outliers/PJ_1outliers.txt \
	-o ../results/biallelic_outliers/rad_locus/pj_1outlier_loci.fasta 
	
#~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
#	-c ../results/stacks/batch_1.catalog.tags.tsv \
#	-l ../results/biallelic_outliers/viability_outliers.txt \
#	-o ../results/biallelic_outliers/rad_locus/viability_outlier_loci.fasta 
	
~/sf_ubuntushare/fasta_from_stacks_catalog/fasta_from_stacks_catalog/fasta_from_catalog \
	-c ../results/stacks/batch_1.catalog.tags.tsv \
	-l ../results/biallelic_outliers/shared1.txt \
	-o ../results/biallelic_outliers/rad_locus/shared_1outlier_loci.fasta 
