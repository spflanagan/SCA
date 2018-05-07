#!/bin/bash

BIALLELIC=true
HAPLOTYPES=false
SIMULATION=true

LOCI="50 100 150 200 300 400 800 1600"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


cd $DIR

if [ "$HAPLOTYPES" = true ]; then
	PROGDIR="../results/parentage_haplotypes"
	cd $PROGDIR

	for L in ${LOCI}
	do
		for r in {1..10}
		do
			m=$((${L}/2))
			sed -e "s/dradPrunedHaps50_1/dradPrunedHaps${L}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" maternity_template.crv > dradPrunedHaps${L}_${r}_mat.crv
			sed -e "s/dradPrunedHaps50_1/dradPrunedHaps${L}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" paternity_template.crv > dradPrunedHaps${L}_${r}_pat.crv
		done
	done
fi

cd $DIR

if [ "$BIALLELIC" = true ]; then
	PROGDIR="../results/parentage_biallelic"
	cd $PROGDIR

	for L in ${LOCI}
	do
		for r in {1..10}
		do
			m=$((${L}/2))
			sed -e "s/dradPruned50_1/dradPruned${L}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" maternity_template.crv > dradPruned${L}_${r}_mat.crv
			sed -e "s/dradPruned50_1/dradPruned${L}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" paternity_template.crv > dradPruned${L}_${r}_pat.crv
		done
	done
fi

cd $DIR

if [ "$SIMULATION" = true ]; then

	PROGDIR="../results/parentsim"
	cd $PROGDIR
	NFEM="50 100 150"
	NSNP="1 4"
	for r in {1..10}
	do
		for L in ${LOCI}
		do
			for f in ${NFEM}
			do
				for s in ${NSNP}
				do
					m=$((${L}/2))
					sed -e "s/parentsim_L50S1F50_1/parentsim_L${L}S${s}F${f}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" -e "s/NCandidateFemales=50/NCandidateFemales=${f}/g" parentsim_maternity_template.crv > parentsim_L${L}S${s}F${f}_${r}_mat.crv
					sed -e "s/parentsim_L50S1F50_1/parentsim_L${L}S${s}F${f}_${r}/g" -e "s/NLoci=50/NLoci=${L}/g" -e "s/MinTypedLoci=25/MinTypedLoci=${m}/g" -e "s/NCandidateMales=50/NCandidateMales=${f}/g" parentsim_paternity_template.crv > parentsim_L${L}S${s}F${f}_${r}_pat.crv
				done
			done
		done
	done
fi

