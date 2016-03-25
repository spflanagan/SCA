#!/bin/bash
cd ../results/biallelic

bowtie2 --sensitive -x ~/sf_ubuntushare/scovelli_genome/annotated_scovelli \
	-U top1.out.radloc.fasta -f \
	-S top1.out.radloc_align.sam -t -p 6
