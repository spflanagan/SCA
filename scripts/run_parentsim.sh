#!/bin/bash

# Author: Sarah P. Flanagan
# Date: 12 March 2018
# This script runs multiple reps of the parentage simulation program

###---SET THESE PARAMETERS---###
NUMREPS=10
SNPS="10 50 100 150 200 300 400 800 1600"
NFEM="50 100 500 2500"

DATE=`date +%Y%m%d`

###---RUN FROM WHEREVER---###
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../programs/parentage_sim"

cd $DIR
cd $PROGDIR

### --- RUN THE PARAMETER COMBINATIONS --- ###
echo "Running ${NUMREPS} reps of the parentage simulation program"
echo "Testing parameters ${SNPS} and ${NFEM}"
echo "The program will run in the background."
echo "Check the status with htop or by looking at parentsim_${DATE}.log"


for i in `seq 1 $NUMREPS`
do
	for j in ${SNPS}
	do
		for jj in ${NFEM}
		do
			./parentsim -S ${j} -F ${jj} -M ${jj} -o ../../results/parentsim/parentsim_S${j}F${jj}_${i}
		done
	done
done >> ../../parentsim_${DATE}.log 2>&1 &