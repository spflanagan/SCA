#!/bin/bash

cd ~/sf_ubuntushare/SCA/

#run the pre-gatk steps, including bwa alignment
##NOTE: This takes about a week
sh ./scripts/pre_gatk.sh

#run gatk
sh ./scripts/run_gatk.sh

#Now filter vcf output
vcftools --vcf ./results/gatk/genotype_output.vcf --min-alleles 2 --max-alleles 2 --maf 0.05 --max-maf 0.95 --min-meanDP 1.0 --max-meanDP 100.0 --recode --recode-INFO-all --out ./results/gatk/genotype_output_pruned

cp ./results/gatk/genotype_output_pruned.recode.vcf ./results/monnahan/
cd ./results/monnahan/
#Check depth
python het_v_depth.mod.py gatkgenotype_output_pruned.recode.vcf

