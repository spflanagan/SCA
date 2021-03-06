#!/bin/bash/
cd ../samples/

files="sample_FEM082
sample_PRM177
sample_PRM177-1"

for file in $files
do
	bowtie2 --sensitive -x scovelli_allpaths \
	-S ${file}_align.sam -t \
	-U ${file}.fq \
	-p 6
done
