#!/bin/bash
cd ~/sf_ubuntushare/ovary_analysis/results/biallelic

vcftools --vcf ../stacks/batch_1.vcf --keep ../extract_from_vcf.txt --extract-FORMAT-info GT --out fem
vcftools --vcf biallelic_maternal.vcf --extract-FORMAT-info GT --out biallelic_maternal
