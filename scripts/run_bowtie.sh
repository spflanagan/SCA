#!/bin/bash

cd ../samples/
files="sample_FEM002
sample_FEM004
sample_FEM010
sample_FEM001
sample_FEM005
sample_FEM006
sample_FEM008
sample_FEM009
sample_FEM011
sample_FEM013
sample_FEM014
sample_FEM015
sample_FEM016
sample_FEM017
sample_FEM018
sample_FEM019
sample_FEM020
sample_FEM021
sample_FEM022
sample_FEM023
sample_FEM024
sample_FEM025
sample_FEM026
sample_FEM027
sample_FEM028
sample_FEM029
sample_FEM030
sample_FEM031
sample_FEM032
sample_FEM033
sample_FEM034
sample_FEM035
sample_FEM036
sample_FEM037
sample_FEM039
sample_FEM041
sample_FEM046
sample_FEM047
sample_FEM051
sample_FEM052
sample_FEM054
sample_FEM056
sample_FEM058
sample_FEM064
sample_FEM066
sample_FEM067
sample_FEM072
sample_FEM074
sample_FEM076
sample_FEM077
sample_FEM078
sample_FEM080
sample_FEM081
sample_FEM083
sample_FEM084
sample_FEM087
sample_NPM005
sample_NPM006
sample_NPM007
sample_NPM008
sample_NPM010
sample_NPM1128
sample_NPM011
sample_NPM012
sample_OFF001
sample_OFF004
sample_OFF005
sample_OFF006
sample_OFF007
sample_OFF008
sample_OFF009
sample_OFF010
sample_OFF011
sample_OFF012
sample_OFF013
sample_OFF014
sample_OFF015
sample_OFF016
sample_ROFF016
sample_OFF017
sample_OFF018
sample_OFF020
sample_OFF022
sample_OFF024
sample_OFF025
sample_OFF026
sample_OFF027
sample_ROFF027
sample_OFF028
sample_OFF029
sample_OFF030
sample_OFF031
sample_OFF032
sample_ROFF032
sample_OFF033
sample_OFF034
sample_OFF035
sample_OFF036
sample_OFF037
sample_OFF038
sample_OFF039
sample_OFF041
sample_OFF042
sample_OFF043
sample_OFF044
sample_OFF045
sample_OFF046
sample_OFF047
sample_OFF049
sample_OFF050
sample_OFF051
sample_OFF052
sample_OFF053
sample_OFF054
sample_OFF055
sample_OFF056
sample_OFF057
sample_OFF058
sample_OFF059
sample_OFF060
sample_OFF061
sample_OFF064
sample_OFF066
sample_OFF067
sample_OFF068
sample_OFF070
sample_OFF071
sample_OFF072
sample_OFF073
sample_OFF074
sample_OFF075
sample_OFF076
sample_OFF077
sample_OFF078
sample_OFF079
sample_OFF080
sample_OFF081
sample_OFF083
sample_OFF084
sample_OFF085
sample_OFF08623
sample_OFF086
sample_OFF088
sample_OFF089
sample_OFF090
sample_OFF091
sample_OFF092
sample_OFF093
sample_OFF094
sample_OFF095
sample_OFF096
sample_OFF097
sample_OFF100
sample_OFF101
sample_OFF102
sample_OFF103
sample_OFF105
sample_OFF106
sample_OFF110
sample_OFF111
sample_OFF112
sample_OFF113
sample_OFF114
sample_OFF115
sample_OFF117
sample_OFF118
sample_OFF119
sample_OFF120
sample_OFF121
sample_OFF122
sample_OFF123
sample_OFF124
sample_OFF125
sample_OFF126
sample_OFF127
sample_OFF134
sample_OFF135
sample_OFF136
sample_OFF137
sample_OFF138
sample_OFF139
sample_OFF140
sample_OFF142
sample_OFF143
sample_OFF144
sample_OFF145
sample_OFF146
sample_OFF149
sample_OFF150
sample_OFF151
sample_OFF152
sample_OFF153
sample_OFF154
sample_OFF155
sample_OFF156
sample_OFF157
sample_OFF158
sample_OFF159
sample_OFF160
sample_OFF161
sample_OFF165
sample_OFF166
sample_OFF167
sample_OFF168
sample_OFF169
sample_OFF170
sample_OFF171
sample_OFF172
sample_OFF173
sample_OFF174
sample_OFF175
sample_OFF176
sample_OFF177
sample_OFF178
sample_OFF179
sample_OFF181
sample_OFF182
sample_OFF183
sample_OFF184
sample_OFF185
sample_OFF186
sample_OFF187
sample_OFF188
sample_OFF189
sample_PRM001
sample_PRM002
sample_PRM003
sample_PRM005
sample_PRM006
sample_PRM007
sample_PRM009
sample_PRM010
sample_PRM011
sample_PRM012
sample_PRM013
sample_PRM014
sample_PRM015
sample_PRM016
sample_PRM017
sample_PRM018
sample_PRM019
sample_PRM022
sample_PRM023
sample_PRM025
sample_PRM026
sample_PRM028
sample_PRM041
sample_PRM042
sample_PRM043
sample_PRM044
sample_PRM045
sample_PRM046
sample_PRM047
sample_PRM048
sample_PRM049
sample_PRM050
sample_PRM051
sample_PRM053
sample_PRM054
sample_PRM055
sample_PRM056
sample_PRM057
sample_PRM058
sample_PRM059
sample_PRM060
sample_PRM061
sample_PRM062
sample_PRM064
sample_PRM065
sample_PRM066
sample_PRM067
sample_PRM068
sample_PRM069
sample_PRM070
sample_PRM071
sample_PRM072
sample_PRM073
sample_PRM074
sample_PRM076
sample_PRM077
sample_PRM078
sample_PRM079
sample_PRM080
sample_PRM081
sample_PRM082
sample_PRM083
sample_PRM084
sample_PRM085
sample_PRM086-23
sample_PRM086R
sample_PRM087
sample_PRM088
sample_PRM089
sample_PRM090
sample_PRM091
sample_PRM092
sample_PRM093
sample_PRM096
sample_PRM097
sample_PRM098
sample_PRM099
sample_PRM100
sample_PRM101
sample_PRM102
sample_PRM103
sample_PRM104
sample_PRM105
sample_PRM106
sample_PRM107
sample_PRM108
sample_PRM109
sample_PRM110
sample_PRM111
sample_PRM112
sample_PRM113
sample_PRM114
sample_PRM115
sample_PRM117
sample_PRM118
sample_PRM119
sample_PRM120
sample_PRM121
sample_PRM122
sample_PRM124
sample_PRM125
sample_PRM126
sample_PRM127
sample_PRM128
sample_PRM129
sample_PRM130
sample_PRM131
sample_PRM132
sample_PRM133
sample_PRM134
sample_PRM136
sample_PRM137
sample_PRM139
sample_PRM140
sample_PRM143
sample_PRM144
sample_PRM145
sample_PRM146
sample_PRM147
sample_PRM148
sample_PRM149
sample_PRM150
sample_PRM151
sample_PRM152
sample_PRM153
sample_PRM154
sample_PRM155
sample_PRM156
sample_PRM159
sample_PRM160
sample_PRM161
sample_PRM162
sample_PRM163
sample_PRM164
sample_PRM165
sample_PRM166
sample_PRM167
sample_PRM168
sample_PRM169
sample_PRM171
sample_PRM172
sample_PRM173
sample_PRM174
sample_PRM175
sample_PRM176
sample_PRM177
sample_PRM177
sample_PRM178
sample_PRM179
sample_PRM181
sample_PRM182
sample_PRM183
sample_PRM184
sample_PRM185
sample_PRM186
sample_PRM187
sample_PRM188
sample_PRM189"

for file in $files
do
	bowtie2 --sensitive -x scovelli_allpaths \
	-S ${file}_align.sam -t \
	-U ${file}.fq \
	-p 6
done
