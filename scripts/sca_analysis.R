#Author: Sarah P. Flanagan
#Date: 10 August 2015
#Purpose: Analyze SCA data from Stacks

rm(list=ls())

#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################
sca.map<-read.table("E://ubuntushare//SCA//sca_popmap.txt")
sex<-as.character(sca.map[,2])
sex[sex=="FEM"]<-2
sex[sex == "NPM"]<-1
sex[sex == "PRM"]<-1
sex[sex == "OFF"]<-0
sca.ind.info<-data.frame(
	IndID=sub('sample_(\\w{3}\\d{3})_align','\\1',sca.map[,1]), 
	year=rep(2011),locality=rep(1),
	population=sca.map[,2],sex=sex,age=rep(-9),stay=rep(-9),
	recruits=rep(-9),survival=rep(-9))
write.table(sca.ind.info, "E://ubuntushare//SCA//sca.ind.info.txt", 
	col.names=T,row.names=T,eol='\n',sep='\t', quote=F)
