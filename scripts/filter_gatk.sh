#!/bin/bash

cd results/gatk/

java -jar ~/Programs/GATK/GenomeAnalysisTK.jar \
   -T SelectVariants \
   -R allpaths_cms1.scaff.fa \
   -V genotype_output.vcf \
   -o filtered_bi.vcf \
   -selectType SNP \
   -restrictAllelesTo BIALLELIC 
   
java -jar ~/Programs/GATK/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R allpaths_cms1.scaff.fa \
	-o filtered_output.vcf \
	--variant filtered_bi.vcf \
	--filterName "allelefreq" --filterExpression "AF < 0.05 && AF > 0.95" \
	--filterName "depth" --filterExpression "DP < 1.0 && DP > 100.0" \
	--filterName "allelenum" --filterExpression "AN != 2.0"\
	--setFilteredGtToNocall 

