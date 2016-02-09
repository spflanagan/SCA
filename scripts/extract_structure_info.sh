#!/bin/bash

cd ~/sf_ubuntushare/SCA
FILES="./results/structure/sca/admixture/Results/*_f"

for f in $FILES
do
	echo "Parsing ${f}"
	../scovelli_popgen/saltwater/programs/parse_structure_output/parse_structure_output/parse_structure_output -i ${f}	
done
