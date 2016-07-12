#Author: Sarah P. Flanagan
#Last Updated: 10 June 2016
#Date: 14 February 2016
#Purpose: Analyze relationships between female mating and reproductive 
#success and female morhpometrics.

##################################FUNCTIONS#################################
standardize.trait<-function(vec){
	std.vec<-(vec-mean(vec))/sd(vec)
	return(std.vec)
}

relative.fit<-function(vec){
	rel.vec<-vec/mean(vec)
	return(rel.vec)
}

ci95<-function(vec){
	ci<-c(mean(vec)-(1.96*(sd(vec)/sqrt(length(vec)))),
		mean(vec)+(1.96*(sd(vec)/sqrt(length(vec)))))
	return(ci)
}
####################################ANALYSIS##################################

setwd("B:/ubuntushare/SCA/results/morphometrics/")
fem.morph.dat<-read.csv("scovelli11_morph.csv")
fem.morph.dat$Fish.ID<-gsub("S11F(\\d+)","FEM\\1",fem.morph.dat$Fish.ID)
mal.morph.dat<-read.csv("MaleScovelliPhenotype.csv")
mal.morph.dat$Fish.ID<-gsub("S11PM(\\d+)","PRM\\1",mal.morph.dat$Fish.ID)
mal.morph.dat$Fish.ID<-gsub("S11(NPM\\d+)","\\1",mal.morph.dat$Fish.ID)
mat.dat<-read.delim("../parentage/batemanator_input.rerun.txt")
mat.dat$Fish.ID<-gsub("S11PM(\\d+)","PRM\\1",mat.dat$Fish.ID)
mat.dat$Fish.ID<-gsub("S11(NPM\\d+)","\\1",mat.dat$Fish.ID)



frs.dat<-merge(fem.morph.dat, mat.dat, by="Fish.ID",all.x=T)
frs.dat$TailLength<-frs.dat$StdLength-frs.dat$SVL
frs.dat$HeadLength<-frs.dat$HeadLength-frs.dat$SnoutLength
write.csv(frs.dat,"FemaleRSandMorph.csv")

fem.dat<-frs.dat
fem.dat$NumMates[is.na(fem.dat$NumMates)]<-0
fem.dat$No.Offspring[is.na(fem.dat$No.Offspring)]<-0

ms.num.aov<-lm(NumMates~BandNum,dat=fem.dat)
anova(ms.num.aov)
ms.area.aov<-lm(NumMates~MeanBandArea,dat=fem.dat)
anova(ms.area.aov)
rs.num.aov<-lm(No.Offspring~BandNum,dat=fem.dat)
anova(rs.num.aov)
ms.area.aov<-lm(No.Offspring~MeanBandArea,dat=fem.dat)
anova(rs.area.aov)
#none of these are significant.

#################SELECTION DIFFERENTIALS#############################
frs.dat<-frs.dat[!is.na(frs.dat$NumMates),]
std.fem<-data.frame(frs.dat$Fish.ID,apply(
	frs.dat[,c("NumMates","No.Offspring")],2,relative.fit),
	apply(frs.dat[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"HeadLength","BandNum","MeanBandArea")],2,standardize.trait))

mprime<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"HeadLength","BandNum","MeanBandArea")],2,cov, std.fem$NumMates)
sprime<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"HeadLength","BandNum","MeanBandArea")],2,cov, std.fem$No.Offspring)
#sprime2<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
#		"HeadLength","BandNum","MeanBandArea")],2,cov, std.fem$SurvivingOff)
sel.diff<-rbind(mprime,sprime)
rownames(sel.diff)<-c("Mating Success", "Reproductive Success")
write.csv(sel.diff, "SelectionDifferentials.csv")
