#!/bin/bash

##### gbstools #####
##**NEED TO UPDATE print() COMMANDS THROUGHTOUT GBSTOOLS**##
#get the restriction sites
#first download withrefm and proto files, and run 'sudo rebaseextract`
restrict -sequence ~/sf_ubuntushare/scovelli_genome/SSC_integrated.fa -enzymes PstI,MboI -outfile ~/sf_ubuntushare/SCA/results/ssc.digest -sitelen 2
cat ssc.digest | python ~/Programs/gbstools/bin/digest_to_bed.py > ssc.digest.bed

#convert restriction site BED to GBS-BED with tabix
#I added MboI to make_gbsbed.py before installing 
cat ssc.digest.bed | make_gbsbed.py > ssc.digest.gbsbed
cat ssc.digest.gbsbed | bgzip -c > ssc.digest.gbsbed.gz
tabix -p bed ssc.digest.gbsbed.gz

#I have a vcf file and no normfactors, so I'll just let them default to 1.0
#apparently had to download scipy
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i both.sub.vcf -o both.lrt.vcf

#Now with drad
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i drad.sub.vcf -o drad.lrt.vcf

#and orad
python2.7 ~/Programs/gbstools/bin/polymorphism_test.py -i orad.sub.vcf -o orad.lrt.vcf