#!/bin/bash

##### gbstools #####
##**NEED TO UPDATE print() COMMANDS THROUGHTOUT GBSTOOLS**##
#get the restriction sites
#first download withrefm and proto files, and run 'sudo rebaseextract`
#restrict -sequence ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -enzymes PstI,MboI -outfile ~/sf_ubuntushare/SCA/results/ssc.ddigest -sitelen 2
#cat ssc.ddigest | python2.7 ~/Programs/gbstools/bin/digest_to_bed.py > ssc.ddigest.bed
mv ssc.digest ssc.ddigest
mv ssc.ddigest.bed ssc.ddigest.bed
restrict -sequence ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -enzymes PstI -outfile ~/sf_ubuntushare/SCA/results/ssc.sdigest -sitelen 2
cat ssc.sdigest | python2.7 ~/Programs/gbstools/bin/digest_to_bed.py > ssc.sdigest.bed
#convert restriction site BED to GBS-BED with tabix
#I added MboI to make_gbsbed.py before installing 
cat ssc.ddigest.bed | python 2.7 ~/Programs/gbstools/bin/make_gbsbed.py > ssc.ddigest.gbsbed
cat ssc.ddigest.gbsbed | bgzip -c > ssc.ddigest.gbsbed.gz
tabix -p bed ssc.ddigest.gbsbed.gz
cat ssc.sdigest.bed | python 2.7 ~/Programs/gbstools/bin/make_gbsbed.py > ssc.sdigest.gbsbed
cat ssc.sdigest.gbsbed | bgzip -c > ssc.sdigest.gbsbed.gz
tabix -p bed ssc.sdigest.gbsbed.gz

#I have a vcf file and no normfactors, so I'll just let them default to 1.0
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i both.sub.vcf -o both.lrt.vcf

#Now with drad
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i drad.sub.vcf -o drad.lrt.vcf

#and orad
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i orad.sub.vcf -o orad.lrt.vcf