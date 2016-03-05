#!/bin/bash

cd results/gatk/
java -jar ~/Programs/GATK/GenomeAnalysisTK.jar \
   -T VariantFiltration \
   -R allpaths_cms1.scaff.fa \
   -o filtered_output.vcf \
   --variant genotype_output.vcf \
   --filterName allelefreq --filterExpression "AF < 0.05 && AF < 0.95" \
   --filterName depth --filterExpression "DP >1 && DP <= 100" \
   --filterName allelenum --filterExpression "AN == 2"

