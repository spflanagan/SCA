#!/bin/bash

cd ../samples/PstI/
files="orad_FEM007
orad_FEM054
orad_FEM071
orad_PRM033
orad_PRM095
orad_FEM012
orad_FEM055
orad_FEM073
orad_PRM034
orad_PRM116
orad_FEM038
orad_FEM057
orad_FEM075
orad_PRM035
orad_PRM123
orad_FEM040
orad_FEM059
orad_FEM079
orad_PRM036
orad_PRM135
orad_FEM042
orad_FEM060
orad_FEM085
orad_PRM037
orad_PRM138
orad_FEM043
orad_FEM061
orad_FEM086
orad_PRM038
orad_PRM141
orad_FEM044
orad_FEM062
orad_PRM024
orad_PRM039
orad_PRM142
orad_FEM045
orad_FEM063
orad_PRM027
orad_PRM040
orad_PRM157
orad_FEM048
orad_FEM065
orad_PRM029
orad_PRM052
orad_PRM158
orad_FEM049
orad_FEM068
orad_PRM030
orad_PRM063
orad_PRM170
orad_FEM050
orad_FEM069
orad_PRM031
orad_PRM075
orad_PRM180
orad_FEM053
orad_FEM070
orad_PRM032
orad_PRM094"

for file in $files
do
	bowtie2 --sensitive -x ../../../scovelli_genome/ssc_chromonome \
	-S ${file}.sam -t \
	-U ${file}.fq \
	-p 4
done
