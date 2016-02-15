#Author: Sarah P. Flanagan
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
morph.dat<-read.csv("scovelli11_morph.csv")
morph.dat$Fish.ID<-gsub("S11F(\\d+)","FEM\\1",morph.dat$Fish.ID)
mat.dat<-read.csv("../parentage/gen1600_8.maternity.sig.csv")
mrs.dat<-read.csv("maleRS.csv")
mrs.dat$Fish.ID<-gsub("S11PM(\\d+)","PRM\\1",mrs.dat$Fish.ID)

frs.dat<-merge(mrs.dat, mat.dat, by.x="Fish.ID",by.y="Father.ID",all.y=T)
frs.dat[frs.dat$Fish.ID=="PRM035-1",c("SurvivingOff","Reduced")]<-
	mrs.dat[mrs.dat$Fish.ID=="PRM035",2:3]
frs.dat[frs.dat$Fish.ID=="PRM177-1",c("SurvivingOff","Reduced")]<-
	mrs.dat[mrs.dat$Fish.ID=="PRM177",2:3]

frs.dat<-frs.dat[,c("Fish.ID","SurvivingOff","Reduced", "Candidate.mother.ID")]
fem.dat<-data.frame(NumMates=summary(frs.dat$Candidate.mother.ID),
	SurvivingOff=tapply(frs.dat$SurvivingOff,frs.dat$Candidate.mother.ID,sum),
	ReducedEmbryos=tapply(frs.dat$Reduced, frs.dat$Candidate.mother.ID,sum))
fem.dat$FemID<-rownames(fem.dat)
fem.dat<-merge(fem.dat, morph.dat,by.x="FemID",by.y="Fish.ID",all.x=T)
fem.dat[fem.dat$FemID=="FEM001",5:11]<-
	morph.dat[morph.dat$Fish.ID=="FEM001D",2:8]
fem.dat[fem.dat$FemID=="FEM054-1",5:11]<-
	morph.dat[morph.dat$Fish.ID=="FEM054",2:8]
write.csv(fem.dat,"FemaleRSandMorph.csv")

fem.dat$TailLength<-fem.dat$StdLength-fem.dat$SVL
fem.dat$HeadLength<-fem.dat$HeadLength-fem.dat$SnoutLength
fem.dat$TotalOff<-fem.dat$SurvivingOff+fem.dat$ReducedEmbryos
std.fem<-data.frame(fem.dat$FemID,apply(
	fem.dat[,c("NumMates","SurvivingOff","TotalOff")],2,relative.fit),
	apply(fem.dat[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"BandNum","MeanBandArea")],2,standardize.trait))

mprime<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"BandNum","MeanBandArea")],2,cov, std.fem$NumMates)
sprime<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"BandNum","MeanBandArea")],2,cov, std.fem$TotalOff)
sprime2<-apply(std.fem[,c("SVL","SnoutLength","SnoutDepth","TailLength",
		"BandNum","MeanBandArea")],2,cov, std.fem$SurvivingOff)
sel.diff<-rbind(mprime,sprime,sprime2)
rownames(sel.diff)<-c("Mating Success", "Total Num Embryos", "Surviving Embryos")
write.csv(sel.diff, "SelectionDifferentials.csv")