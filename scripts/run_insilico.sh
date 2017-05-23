#!/bin/bash
cd ../results/insilico/
echo "PANEL 1: PCR bias"
../../programs/insilico_radseq/insilico_radseq -o null -p 0 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr1.notskewed -p 0.01 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr2.notskewed -p 0.02 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr3.notskewed -p 0.03 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr4.notskewed -p 0.04 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr5.notskewed -p 0.05 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o null.skewedAFS -p 0 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr1.skewed -p 0.01 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr2.skewed -p 0.02 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr3.skewed -p 0.03 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr4.skewed -p 0.04 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o pcr5.skewed -p 0.05 -c 1 -a 200 -d 200 -a 1 -b 0
echo "PANEL 2: Polymorphism"
../../programs/insilico_radseq/insilico_radseq -o u10E8Ne5000 -p 0 -mu 10E8 -n 5000 -c 0.1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o u10E8Ne10000 -p 0 -mu 10E8 -n 10000 -c 0.1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o u10E8Ne20000 -p 0 -mu 10E8 -n 20000 -c 0.1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o u10E7Ne5000 -p 0 -mu 10E7 -n 5000 -c 0.1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o u10E7Ne10000 -p 0 -mu 10E7 -n 10000 -c 0.1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o u10E7Ne20000 -p 0 -mu 10E7 -n 20000 -c 0.1 -a 200 -d 200 -a 0 -b 0
echo "PANEL 3: Skew, Shearing Bias, and sampling"
../../programs/insilico_radseq/insilico_radseq -o null -p 0 -c 1 -a 200 -d 200 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o shearingbias.notskewed -p 0 -c 1 -a 200 -d 200 -a 0 -b 1
../../programs/insilico_radseq/insilico_radseq -o shearingbias.skewed -p 0 -c 1 -a 200 -d 200 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o null.skewedAFS -p 0 -c 1 -a 200 -d 200 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o null.skewedAFS.asymmetric -p 0 -c 1 -a 60 -d 340 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o null.asymmetric -p 0 -c 1 -a 60 -d 340 -a 0 -b 0
../../programs/insilico_radseq/insilico_radseq -o shearingbias.asymmetric.notskewed -p 0 -c 1 -a 60 -d 340 -a 0 -b 1
../../programs/insilico_radseq/insilico_radseq -o shearingbias.asymmetric.skewed -p 0 -c 1 -a 60 -d 340 -a 1 -b 1
echo "PANEL 4: Combinations"
../../programs/insilico_radseq/insilico_radseq -o comb.pcr1.u10E8Ne10000.asymmetric.skewed.shearingbias -p 0.01 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr2.u10E8Ne10000.asymmetric.skewed.shearingbias -p 0.02 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr3.u10E8Ne10000.asymmetric.skewed.shearingbias -p 0.03 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr1.u10E7Ne10000.asymmetric.skewed.shearingbias -p 0.01 -mu 10E7 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr2.u10E7Ne10000.asymmetric.skewed.shearingbias -p 0.02 -mu 10E7 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr3.u10E7Ne10000.asymmetric.skewed.shearingbias -p 0.03 -mu 10E7 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr1.u10E8Ne10000.symmetric.skewed.shearingbias -p 0.01 -mu 10E8 -n 10000 -c 0.1 -a 200 -d 200 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr2.u10E8Ne10000.symmetric.skewed.shearingbias -p 0.02 -mu 10E8 -n 10000 -c 0.1 -a 200 -d 200 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr3.u10E8Ne10000.symmetric.skewed.shearingbias -p 0.03 -mu 10E8 -n 10000 -c 0.1 -a 200 -d 200 -a 1 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr1.u10E8Ne10000.asymmetric.shearingbias -p 0.01 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 0 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr2.u10E8Ne10000.asymmetric.shearingbias -p 0.02 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 0 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr3.u10E8Ne10000.asymmetric.shearingbias -p 0.03 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 0 -b 1
../../programs/insilico_radseq/insilico_radseq -o comb.pcr1.u10E8Ne10000.asymmetric.skewed -p 0.01 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o comb.pcr2.u10E8Ne10000.asymmetric.skewed -p 0.02 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 0
../../programs/insilico_radseq/insilico_radseq -o comb.pcr3.u10E8Ne10000.asymmetric.skewed -p 0.03 -mu 10E8 -n 10000 -c 0.1 -a 60 -d 340 -a 1 -b 0