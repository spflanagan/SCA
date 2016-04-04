#!/bin/bash
cd ../results/biallelic_outliers

#bowtie2 --sensitive -x ~/sf_ubuntushare/scovelli_genome/annotated_scovelli \
#	-U rad_locus/aj_outlier_loci.fasta -f \
#	-S align_to_annotated/aj_outlier_loci_align.sam -t -p 6
	
bowtie2 --sensitive -x ~/sf_ubuntushare/scovelli_genome/annotated_scovelli \
	-U rad_locus/fm_1outlier_loci.fasta -f \
	-S align_to_annotated/fm_1outlier_loci_align.sam -t -p 6
	
bowtie2 --sensitive -x ~/sf_ubuntushare/scovelli_genome/annotated_scovelli \
	-U rad_locus/mo_1outlier_loci.fasta -f \
	-S align_to_annotated/mo_1outlier_loci_align.sam -t -p 6
	
bowtie2 --sensitive -x ~/sf_ubuntushare/scovelli_genome/annotated_scovelli \
	-U rad_locus/pj_1outlier_loci.fasta -f \
	-S align_to_annotated/pj_1outlier_loci_align.sam -t -p 6
