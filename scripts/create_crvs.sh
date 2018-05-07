#!/bin/bash

LOCI="50 100 150 200 300 400 800 1600"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROGDIR="../results/parentage_haplotypes"

cd $DIR
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

PROGDIR="../parentage_biallelic"

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