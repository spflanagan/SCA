#!/bin/bash

cd ~/sf_ubuntushare/SCA/results/gatk/
java -jar ~/Programs/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs \
	-R allpaths_cms1.scaff.fa \
	-V FEM001.g.vcf \
	-V FEM002.g.vcf \
	-V FEM004.g.vcf \
	-V FEM005.g.vcf \
	-V FEM006.g.vcf \
	-V FEM008.g.vcf \
	-V FEM009.g.vcf \
	-V FEM010.g.vcf \
	-V FEM011.g.vcf \
	-V FEM013.g.vcf \
	-V FEM014.g.vcf \
	-V FEM015.g.vcf \
	-V FEM016.g.vcf \
	-V FEM017.g.vcf \
	-V FEM018.g.vcf \
	-V FEM019.g.vcf \
	-V FEM020.g.vcf \
	-V FEM021.g.vcf \
	-V FEM022.g.vcf \
	-V FEM023.g.vcf \
	-V FEM024.g.vcf \
	-V FEM025.g.vcf \
	-V FEM026.g.vcf \
	-V FEM027.g.vcf \
	-V FEM028.g.vcf \
	-V FEM029.g.vcf \
	-V FEM030.g.vcf \
	-V FEM031.g.vcf \
	-V FEM032.g.vcf \
	-V FEM033.g.vcf \
	-V FEM034.g.vcf \
	-V FEM035.g.vcf \
	-V FEM036.g.vcf \
	-V FEM037.g.vcf \
	-V FEM039.g.vcf \
	-V FEM041.g.vcf \
	-V FEM046.g.vcf \
	-V FEM047.g.vcf \
	-V FEM051.g.vcf \
	-V FEM052.g.vcf \
	-V FEM054-1.g.vcf \
	-V FEM056.g.vcf \
	-V FEM058.g.vcf \
	-V FEM064.g.vcf \
	-V FEM066.g.vcf \
	-V FEM067.g.vcf \
	-V FEM072.g.vcf \
	-V FEM074.g.vcf \
	-V FEM076.g.vcf \
	-V FEM077.g.vcf \
	-V FEM078.g.vcf \
	-V FEM080.g.vcf \
	-V FEM081.g.vcf \
	-V FEM082.g.vcf \
	-V FEM083.g.vcf \
	-V FEM084.g.vcf \
	-V FEM087.g.vcf \
	-V NPM005.g.vcf \
	-V NPM006.g.vcf \
	-V NPM007.g.vcf \
	-V NPM008.g.vcf \
	-V NPM010.g.vcf \
	-V NPM011.g.vcf \
	-V NPM012.g.vcf \
	-V NPM1128.g.vcf \
	-V OFF001.g.vcf \
	-V OFF004.g.vcf \
	-V OFF005.g.vcf \
	-V OFF006.g.vcf \
	-V OFF007.g.vcf \
	-V OFF008.g.vcf \
	-V OFF009.g.vcf \
	-V OFF010.g.vcf \
	-V OFF011.g.vcf \
	-V OFF012.g.vcf \
	-V OFF013.g.vcf \
	-V OFF014.g.vcf \
	-V OFF015.g.vcf \
	-V OFF016.g.vcf \
	-V OFF017.g.vcf \
	-V OFF018.g.vcf \
	-V OFF020.g.vcf \
	-V OFF022.g.vcf \
	-V OFF024.g.vcf \
	-V OFF025.g.vcf \
	-V OFF026.g.vcf \
	-V OFF027.g.vcf \
	-V OFF028.g.vcf \
	-V OFF029.g.vcf \
	-V OFF030.g.vcf \
	-V OFF031.g.vcf \
	-V OFF032.g.vcf \
	-V OFF033.g.vcf \
	-V OFF034.g.vcf \
	-V OFF035.g.vcf \
	-V OFF036.g.vcf \
	-V OFF037.g.vcf \
	-V OFF038.g.vcf \
	-V OFF039.g.vcf \
	-V OFF041.g.vcf \
	-V OFF042.g.vcf \
	-V OFF043.g.vcf \
	-V OFF044.g.vcf \
	-V OFF045.g.vcf \
	-V OFF046.g.vcf \
	-V OFF047.g.vcf \
	-V OFF049.g.vcf \
	-V OFF050.g.vcf \
	-V OFF051.g.vcf \
	-V OFF052.g.vcf \
	-V OFF053.g.vcf \
	-V OFF054.g.vcf \
	-V OFF055.g.vcf \
	-V OFF056.g.vcf \
	-V OFF057.g.vcf \
	-V OFF058.g.vcf \
	-V OFF059.g.vcf \
	-V OFF060.g.vcf \
	-V OFF061.g.vcf \
	-V OFF064.g.vcf \
	-V OFF066.g.vcf \
	-V OFF067.g.vcf \
	-V OFF068.g.vcf \
	-V OFF070.g.vcf \
	-V OFF071.g.vcf \
	-V OFF072.g.vcf \
	-V OFF073.g.vcf \
	-V OFF074.g.vcf \
	-V OFF075.g.vcf \
	-V OFF076.g.vcf \
	-V OFF077.g.vcf \
	-V OFF078.g.vcf \
	-V OFF079.g.vcf \
	-V OFF080.g.vcf \
	-V OFF081.g.vcf \
	-V OFF083.g.vcf \
	-V OFF084.g.vcf \
	-V OFF085.g.vcf \
	-V OFF08623.g.vcf \
	-V OFF086.g.vcf \
	-V OFF088.g.vcf \
	-V OFF089.g.vcf \
	-V OFF090.g.vcf \
	-V OFF091.g.vcf \
	-V OFF092.g.vcf \
	-V OFF093.g.vcf \
	-V OFF094.g.vcf \
	-V OFF095.g.vcf \
	-V OFF096.g.vcf \
	-V OFF097.g.vcf \
	-V OFF100.g.vcf \
	-V OFF101.g.vcf \
	-V OFF102.g.vcf \
	-V OFF103.g.vcf \
	-V OFF105.g.vcf \
	-V OFF106.g.vcf \
	-V OFF110.g.vcf \
	-V OFF111.g.vcf \
	-V OFF112.g.vcf \
	-V OFF113.g.vcf \
	-V OFF114.g.vcf \
	-V OFF115.g.vcf \
	-V OFF117.g.vcf \
	-V OFF118.g.vcf \
	-V OFF119.g.vcf \
	-V OFF120.g.vcf \
	-V OFF121.g.vcf \
	-V OFF122.g.vcf \
	-V OFF123.g.vcf \
	-V OFF124.g.vcf \
	-V OFF125.g.vcf \
	-V OFF126.g.vcf \
	-V OFF127.g.vcf \
	-V OFF134.g.vcf \
	-V OFF135.g.vcf \
	-V OFF136.g.vcf \
	-V OFF137.g.vcf \
	-V OFF138.g.vcf \
	-V OFF139.g.vcf \
	-V OFF140.g.vcf \
	-V OFF142.g.vcf \
	-V OFF143.g.vcf \
	-V OFF144.g.vcf \
	-V OFF145.g.vcf \
	-V OFF146.g.vcf \
	-V OFF149.g.vcf \
	-V OFF150.g.vcf \
	-V OFF151.g.vcf \
	-V OFF152.g.vcf \
	-V OFF153.g.vcf \
	-V OFF154.g.vcf \
	-V OFF155.g.vcf \
	-V OFF156.g.vcf \
	-V OFF157.g.vcf \
	-V OFF158.g.vcf \
	-V OFF159.g.vcf \
	-V OFF160.g.vcf \
	-V OFF161.g.vcf \
	-V OFF165.g.vcf \
	-V OFF166.g.vcf \
	-V OFF167.g.vcf \
	-V OFF168.g.vcf \
	-V OFF169.g.vcf \
	-V OFF170.g.vcf \
	-V OFF171.g.vcf \
	-V OFF172.g.vcf \
	-V OFF173.g.vcf \
	-V OFF174.g.vcf \
	-V OFF175.g.vcf \
	-V OFF176.g.vcf \
	-V OFF177.g.vcf \
	-V OFF178.g.vcf \
	-V OFF179.g.vcf \
	-V OFF181.g.vcf \
	-V OFF182.g.vcf \
	-V OFF183.g.vcf \
	-V OFF184.g.vcf \
	-V OFF185.g.vcf \
	-V OFF186.g.vcf \
	-V OFF187.g.vcf \
	-V OFF188.g.vcf \
	-V OFF189.g.vcf \
	-V PRM001.g.vcf \
	-V PRM002.g.vcf \
	-V PRM003.g.vcf \
	-V PRM005.g.vcf \
	-V PRM006.g.vcf \
	-V PRM007.g.vcf \
	-V PRM009.g.vcf \
	-V PRM010.g.vcf \
	-V PRM011.g.vcf \
	-V PRM012.g.vcf \
	-V PRM013.g.vcf \
	-V PRM014.g.vcf \
	-V PRM015.g.vcf \
	-V PRM016.g.vcf \
	-V PRM017.g.vcf \
	-V PRM018.g.vcf \
	-V PRM019.g.vcf \
	-V PRM022.g.vcf \
	-V PRM023.g.vcf \
	-V PRM025.g.vcf \
	-V PRM026.g.vcf \
	-V PRM028.g.vcf \
	-V PRM041.g.vcf \
	-V PRM042.g.vcf \
	-V PRM043.g.vcf \
	-V PRM044.g.vcf \
	-V PRM045.g.vcf \
	-V PRM046.g.vcf \
	-V PRM047.g.vcf \
	-V PRM048.g.vcf \
	-V PRM049.g.vcf \
	-V PRM050.g.vcf \
	-V PRM051.g.vcf \
	-V PRM053.g.vcf \
	-V PRM054.g.vcf \
	-V PRM055.g.vcf \
	-V PRM056.g.vcf \
	-V PRM057.g.vcf \
	-V PRM058.g.vcf \
	-V PRM059.g.vcf \
	-V PRM060.g.vcf \
	-V PRM061.g.vcf \
	-V PRM062.g.vcf \
	-V PRM064.g.vcf \
	-V PRM065.g.vcf \
	-V PRM066.g.vcf \
	-V PRM067.g.vcf \
	-V PRM068.g.vcf \
	-V PRM069.g.vcf \
	-V PRM070.g.vcf \
	-V PRM071.g.vcf \
	-V PRM072.g.vcf \
	-V PRM073.g.vcf \
	-V PRM074.g.vcf \
	-V PRM076.g.vcf \
	-V PRM077.g.vcf \
	-V PRM078.g.vcf \
	-V PRM079.g.vcf \
	-V PRM080.g.vcf \
	-V PRM081.g.vcf \
	-V PRM082.g.vcf \
	-V PRM083.g.vcf \
	-V PRM084.g.vcf \
	-V PRM085.g.vcf \
	-V PRM086-23.g.vcf \
	-V PRM086R.g.vcf \
	-V PRM087.g.vcf \
	-V PRM088.g.vcf \
	-V PRM089.g.vcf \
	-V PRM090.g.vcf \
	-V PRM091.g.vcf \
	-V PRM092.g.vcf \
	-V PRM093.g.vcf \
	-V PRM096.g.vcf \
	-V PRM097.g.vcf \
	-V PRM098.g.vcf \
	-V PRM099.g.vcf \
	-V PRM100.g.vcf \
	-V PRM101.g.vcf \
	-V PRM102.g.vcf \
	-V PRM103.g.vcf \
	-V PRM104.g.vcf \
	-V PRM105.g.vcf \
	-V PRM106.g.vcf \
	-V PRM107.g.vcf \
	-V PRM108.g.vcf \
	-V PRM109.g.vcf \
	-V PRM110.g.vcf \
	-V PRM111.g.vcf \
	-V PRM112.g.vcf \
	-V PRM113.g.vcf \
	-V PRM114.g.vcf \
	-V PRM115.g.vcf \
	-V PRM117.g.vcf \
	-V PRM118.g.vcf \
	-V PRM119.g.vcf \
	-V PRM120.g.vcf \
	-V PRM121.g.vcf \
	-V PRM122.g.vcf \
	-V PRM124.g.vcf \
	-V PRM125.g.vcf \
	-V PRM126.g.vcf \
	-V PRM127.g.vcf \
	-V PRM128.g.vcf \
	-V PRM129.g.vcf \
	-V PRM130.g.vcf \
	-V PRM131.g.vcf \
	-V PRM132.g.vcf \
	-V PRM133.g.vcf \
	-V PRM134.g.vcf \
	-V PRM135.g.vcf \
	-V PRM136.g.vcf \
	-V PRM137.g.vcf \
	-V PRM139.g.vcf \
	-V PRM140.g.vcf \
	-V PRM143.g.vcf \
	-V PRM144.g.vcf \
	-V PRM145.g.vcf \
	-V PRM146.g.vcf \
	-V PRM147.g.vcf \
	-V PRM148.g.vcf \
	-V PRM149.g.vcf \
	-V PRM150.g.vcf \
	-V PRM151.g.vcf \
	-V PRM152.g.vcf \
	-V PRM153.g.vcf \
	-V PRM154.g.vcf \
	-V PRM155.g.vcf \
	-V PRM156.g.vcf \
	-V PRM159.g.vcf \
	-V PRM160.g.vcf \
	-V PRM161.g.vcf \
	-V PRM162.g.vcf \
	-V PRM163.g.vcf \
	-V PRM164.g.vcf \
	-V PRM165.g.vcf \
	-V PRM166.g.vcf \
	-V PRM167.g.vcf \
	-V PRM168.g.vcf \
	-V PRM169.g.vcf \
	-V PRM171.g.vcf \
	-V PRM172.g.vcf \
	-V PRM173.g.vcf \
	-V PRM174.g.vcf \
	-V PRM175.g.vcf \
	-V PRM176.g.vcf \
	-V PRM177.g.vcf \
	-V PRM178.g.vcf \
	-V PRM179.g.vcf \
	-V PRM181.g.vcf \
	-V PRM182.g.vcf \
	-V PRM183.g.vcf \
	-V PRM184.g.vcf \
	-V PRM185.g.vcf \
	-V PRM186.g.vcf \
	-V PRM187.g.vcf \
	-V PRM188.g.vcf \
	-V PRM189.g.vcf \
	-o genotype_output.vcf
