#!/bin/bash

bowtie2 --sensitive -x ./SCA_samples/scovelli_allpaths -S ./bowtie_gatk/sample_PRM189.bwt.sam -t -U ./SCA_samples/sample_PRM189.fq -p 6 --rg-id sample_PRM189:SM\ssample_PRM189
		
#samtools view -bS ./bowtie_gatk/sample_PRM189.bwt.sam > ./bowtie_gatk/sample_PRM189.bwt.bam
	
	java -jar ~/Programs/picard-tools-1.139/picard.jar CleanSam \
	 INPUT=./bowtie_gatk/sample_PRM189.bwt.bam OUTPUT=./bowtie_gatk/sample_PRM189.bwt_clean.bam
#	java -jar ~/Programs/picard-tools-1.141/picard.jar SortSam \
#	 INPUT=./bowtie_gatk/sample_PRM189.bwt_clean.bam OUTPUT=./bowtie_gatk/sample_PRM189.bwt_sorted.bam \
#	 SORT_ORDER=coordinate
#	java -jar ~/Programs/picard-tools-1.141/picard.jar MarkDuplicates \
#	 INPUT=./bowtie_gatk/sample_PRM189.bwt_sorted.bam OUTPUT=./bowtie_gatk/sample_PRM189.bwt_marked.bam \
#	 METRICS_FILE=./bowtie_gatk/sample_PRM189.bwt.metrics CREATE_INDEX=true
#	java -jar ~/Programs/picard-tools-1.141/picard.jar \
#	 ValidateSamFile INPUT=./bowtie_gatk/sample_PRM189.bwt_marked.bam \
#	 OUTPUT=./bowtie_gatk/sample_PRM189.bwt_validate.log
