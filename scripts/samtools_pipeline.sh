#!/bin/bash

#Title: Re-analyzing RADseq dataset as per reviewer 1's suggestion

##this runs all of the different steps in the samtools alignment
##SUPER thanks to Chris Nice (and Kate Bell) for their scripts, 
##	which helped guide this process.


perl ../scripts/samtools/sam2bam.pl *.sam #ran this on all of the individuals

#run from samples
samtools mpileup -P ILLUMINA -u -g -I -f ../../scovelli_genome/SSC_integrated.fa -b sub60.txt > out1.bcf
bcftools view -N -c -e -g -v -I -d 0.5 -p 0.05 -P full -t 0.001 out1.bcf > ../results/variants50.vcf &
vcftools --vcf ../results/variants50.vcf --out ../results/sta.subset --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --min-meanDP 3 --recode --recode-INFO-all

samtools mpileup -P ILLUMINA -u -g -I -f ../../scovelli_genome/SSC_integrated.fa -b orad60bam.txt > orad.out1.bcf
bcftools view -N -c -e -g -v -I -d 0.5 -p 0.05 -P full -t 0.001 orad.out1.bcf > ../results/orad.variants50.vcf &
vcftools --vcf ../results/orad.variants50.vcf --out ../results/sto.subset --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --min-meanDP 3 --recode --recode-INFO-all

samtools mpileup -P ILLUMINA -u -g -I -f ../../scovelli_genome/SSC_integrated.fa -b drad60bam.txt > drad.out1.bcf
bcftools view -N -c -e -g -v -I -d 0.5 -p 0.05 -P full -t 0.001 drad.out1.bcf > ../results/drad.variants50.vcf &
vcftools --vcf ../results/drad.variants50.vcf --out ../results/std.subset --remove-indels --maf 0.05 --min-alleles 2 --max-alleles 2 --min-meanDP 3 --recode --recode-INFO-all