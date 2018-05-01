#!/bin/bash

# Author: Sarah P. Flanagan
# Date: 12 March 2018
# This script runs multiple reps of the parentage simulation program

###---SET THESE PARAMETERS---###
NUMREPS=10
LOCS="50 100 150 200 300 400 800 1600"
NFEM="50 100 500"
NSNP="1 4"

DATE=`date +%Y%m%d`

###---RUN FROM WHEREVER---###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/parentage_sim"

cd $DIR
cd $PROGDIR

### --- RUN THE PARAMETER COMBINATIONS --- ###
echo "Running ${NUMREPS} reps of the parentage simulation program"
echo "Testing parameters ${LOCS} loci, ${NSNP} SNPs per locus, and ${NFEM} females"
echo "The program will run in the background."
echo "Check the status with htop or by looking at parentsim_${DATE}.log"


for i in `seq 1 $NUMREPS`
do
	for j in ${LOCS}
	do
		for jj in ${NFEM}
		do
			for jjj in ${NSNP}
			do
				./parentsim -L ${j} -S ${jjj} -F ${jj} -M ${jj} -o parentsim_L${j}S${jjj}F${jj}_${i} -d B:\\ubuntushare\\SCA\\results\\parentsim\\ -r ../../results/parentsim/
			done
		done
	done
done >> ../../parentsim_${DATE}.log 2>&1 &
