#!/bin/bash

#Title: Re-analyzing RADseq dataset as per reviewer 1's suggestion

##this runs all of the different steps in the samtools alignment
##SUPER thanks to Chris Nice (and Kate Bell) for their scripts, 
##	which helped guide this process.


perl ../scripts/samtools/sam2bam.pl *.sam #ran this on all of the individuals


samtools mpileup -P ILLUMINA -u -g -I -f ../../scovelli_genome/SSC_integrated.fa -b sub60.txt > ../samples/out1.bcf
bcftools view -N -c -e -g -v -I -d 0.5 -p 0.05 -P full -t 0.001 ../samples/out1.bcf > variants50.vcf &

vcftools --vcf variants50.vcf --out sta.subset --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --min-meanDP 3 --recode --recode-INFO-all