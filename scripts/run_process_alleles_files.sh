#!/bin/bash

list_file="${HOME}/sf_ubuntushare/SCA/fathers_offspring/dad_off_list.txt"

#read file
exec 3<&0
exec 0< $list_file
while read -r line
do	
	#run process_alleles_files
	infile="${HOME}/sf_ubuntushare/SCA/sample_${line}_align.alleles.tsv"
	outfile="${HOME}/sf_ubuntushare/SCA/fathers_offspring/$line.alleles.txt"
	echo "Running process_alleles_files on file $line: in = $infile, out=$outfile"
	./process_alleles_files -i ${infile} -o ${outfile}
done
