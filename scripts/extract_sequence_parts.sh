#!/bin/bash

#run from SCA folder.
./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/shared_2extract.txt -o ./results/biallelic_outliers/rad_region/shared_2extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/mo_1extract.txt -o ./results/biallelic_outliers/rad_region/mo_1extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/fm_1extract.txt -o ./results/biallelic_outliers/rad_region/fm_1extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/aj_extract.txt -o ./results/biallelic_outliers/rad_region/aj_extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/fmmo_extract.txt -o ./results/biallelic_outliers/rad_region/fmmo_extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/ajmo_extract.txt -o ./results/biallelic_outliers/rad_region/ajmo_extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/ajfm_extract.txt -o ./results/biallelic_outliers/rad_region/ajfm_extract.fasta

./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/lrt_fm.txt -o ./results/biallelic_outliers/rad_region/lrt_fm.fasta


./programs/extract_sequence_part/extract_sequence_part -f ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -i ./results/biallelic_outliers/rad_region/lrt_mo.txt -o ./results/biallelic_outliers/rad_region/lrt_mo.fasta