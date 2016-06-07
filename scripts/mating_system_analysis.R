#Author: Sarah P. Flanagan
#Date: 1 June 2016
#Last updated: 1 June 2016
#Purpose: Analyze mating system data from single population of scovelli

rm(list=ls())
#install.packages(“multcomp”, dependencies=TRUE)
library(multcomp)
setwd("E:/ubuntushare/SCA/results/parentage")
dat<-read.delim("batemanator_input.rerun.txt")
fem.dat<-dat[substr(dat$Fish.ID,1,3)=="FEM",]
mal.dat<-dat[substr(dat$Fish.ID,1,3)!="FEM",]
zdat<-read.delim("batemanator_input.rerun_with0s.txt")
femz.dat<-zdat[substr(zdat$Fish.ID,1,3)=="FEM",]

#############################################################################
#########PLOT 
#############################################################################

png("MatingSystem.png",height=10,width=10,res=300,units="in")
#BATEMAN GRADIENTS
par(mar=c(2,2,2,0),oma=c(2,2,0.5,0.5),lwd=1.3)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
#plot males
#plot(mal.dat$NumMates,mal.dat$No.Offspring,pch=17,xlim=c(0,3),ylim=c(-1,70),
#	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
#axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
#axis(2,las=1,pos=0)
#text(x=0.25,y=65,"Males")
#clip(0,1,0,70)
#abline(lm(No.Offspring~NumMates,dat=mal.dat))
#
#mtext("Number of Mates",1,outer=F,line=1,cex=0.85)

#plot females with no zeros
plot(fem.dat$NumMates,fem.dat$No.Offspring,xlim=c(0,3),ylim=c(-1,70),pch=19,
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
#text(x=0.25,y=65,"Females")
clip(1,3,0,70)
#abline(a=0,b=17.69) #use the Bateman gradient slope? 
abline(lm(No.Offspring~NumMates,dat=fem.dat),col="darkorchid4")
clip(0,3,0,70)
abline(lm(No.Offspring~NumMates,dat=femz.dat),lty=2,col="darkorchid4")

mtext("Number of Mates",1,outer=F,line=1,cex=0.85)
mtext("Number of Offspring",2,outer=F,line=1.5,cex=0.85)
#dev.off()
#########PLOT HISTOGRAMS
#png("MS_RS_distributions.png",height=10,width=5,units="in",res=300)
#par(mfrow=c(2,1),mar=c(2,2,2,2),oma=c(2,2,0.5,0.5))
#plot observed female mating success
#hist(mal.dat$NumMates,breaks=seq(-0.5,3.5,1),axes=F,main="",xlab="",ylab="",
#	ylim=c(0,61),col=alpha("dark green",0.5),border=NA)
hist(femz.dat$NumMates,breaks=seq(-0.5,3.5,1),axes=F,main="",xlab="",ylab="",
	ylim=c(0,41),col=alpha("darkorchid4",0.5),border=NA)
axis(1,at=c(-0.5,seq(0,4,1)),lab=c("",seq(0,3,1),""),pos=0,xlim=c(-0.5,4))
axis(2,pos=-0.5,las=1,at=seq(0,40,10),lab=seq(0,40,10))
mtext("Number of Mates",1,outer=F,line=1.5,cex=0.85)
mtext("Number of Individuals",2,outer=F,line=1.5,cex=0.85)
clip(0,4,0,40)
abline(v=mean(femz.dat$NumMates),col="darkorchid4",lty=1)
abline(v=mean(fem.dat$NumMates)-.5,col="darkorchid4",lty=2)



#plot observed female reproductive mating success
hist(mal.dat$No.Offspring, col=alpha("dark green",0.5),border=NA,
	axes=F,main="",xlab="",ylab="",ylim=c(0,61),breaks=seq(0,70,5))
hist(femz.dat$No.Offspring,border=NA,col=alpha("darkorchid4",0.5), add=T,
	breaks=seq(0,70,5))
axis(1,pos=0)
axis(2,pos=0,las=1,at=seq(0,60,10),lab=seq(0,60,10))
mtext("Number of Offspring",1.5,outer=F,line=2,cex=0.85)
legend("topright",
	c("Females","Males","Females With Zeroes","Females Without Zeroes"),
	col=c(alpha("darkorchid4",0.5),alpha("dark green",0.5),
		"darkorchid4","darkorchid4"),
	pch=c(15,15,NA,NA),bty='n',cex=1,lty=c(0,3,2,1))
clip(0,70,0,60)
abline(v=mean(mal.dat$No.Offspring,na.rm=T),col="dark green",lty=3)
abline(v=mean(femz.dat$No.Offspring),col="darkorchid4",lty=1)
abline(v=mean(fem.dat$No.Offspring),col="darkorchid4",lty=2)

dev.off()



#############################################################################
#########BAND SHARING
#############################################################################
all.bs<-read.delim("PolymorphicIn99PercIndsHWE.allcombos.bandsharing.txt")
fo.bs<-read.delim("PolymorphicIn99PercIndsHWE.bandsharing.txt")
all.bs$combo<-paste(all.bs$Father,all.bs$Offspring,sep=".")
all.bs$combo2<-paste(all.bs$Offspring,all.bs$Father,sep=".")

#father-offspring combos
fo.bs$combo<-paste(fo.bs$Father,fo.bs$Offspring,sep=".")

#keep unrelated ones
unr.bs<-all.bs[!(all.bs$combo %in% fo.bs$combo) & 
	!(all.bs$combo2 %in% fo.bs$combo),]
unr.bs$set<-rep("unrelated",nrow(unr.bs))
fo.bs$set<-rep("father-offspring",nrow(fo.bs))

#get combos from maternity analysis
mothers<-read.csv("../parentage/parentage_summary.csv",
	header=T,row.names=NULL)
mothers$combo<-paste(mothers$Offspring.ID,mothers$Candidate.mother.ID,sep=".")
mothers$combo2<-paste(mothers$Candidate.mother.ID,mothers$Offspring.ID,sep=".")

#merge
all.bs<-rbind(unr.bs[,colnames(unr.bs)[colnames(unr.bs)!="combo2"]],fo.bs)
all.bs[all.bs$set=="father-offspring" & all.bs$Shared < 0.2,
	"set"]<-"unrelated"
all.bs[all.bs$combo %in% mothers$combo | all.bs$combo %in% mothers$combo2,
	"set"]<-"mother-offspring"

shared.dat<-data.frame(prop=all.bs$Shared,relationship=all.bs$set,
	type="shared")
incomp.dat<-data.frame(prop=all.bs$Incompatible,relationship=all.bs$set,
	type="incompatible")
all.dat<-rbind(shared.dat,incomp.dat)


#####PLOT
png("band_sharing.png",height=7,width=7,units="in",res=300)
boxplot(all.dat$prop~all.dat$relationship*all.dat$type,axes=F,ylim=c(-0.01,1),
	border=c("dark green","darkorchid4","dodgerblue3"),pch=c(15,17,1))
axis(1,at=c(0,2,5,7),c("","Shared","Incompatible",""),pos=0)
axis(2,pos=0.35,las=1,ylim=c(-0.01,1))
mtext("Proportion of Loci",2,line=2)
legend(x=0.35,y=0.18,
	c("Father-Offspring","Mother-Offspring","Putatively Unrelated"),
	col=c("dark green","darkorchid4","dodgerblue3"),pch=c(15,17,1),bty='n')
dev.off()
#require(ggplot2)
#ggplot(data = all.dat, aes(x=relationship, y=prop)) + 
#	geom_boxplot(aes(fill=type))

#####STATS
sem<-function(x){
	sem<-sd(x)/sqrt(length(x))
}
shared.aov<-lm(prop~relationship,data=shared.dat)
anova(shared.aov)
#Response: prop
#                Df Sum Sq  Mean Sq F value   Pr(>F)    
#relationship     2   0.18 0.092474   20.93 8.18e-10 ***
#Residuals    73533 324.88 0.004418  
shared.tukey<- glht(shared.aov, linfct = mcp(relationship="Tukey"))
summary(shared.tukey)
#Linear Hypotheses:
#                                          Estimate Std. Error t value Pr(>|t|)
#mother-offspring - father-offspring == 0  0.006504   0.014058   0.463  0.88122
#unrelated - father-offspring == 0        -0.033152   0.005835  -5.682  < 0.001
#unrelated - mother-offspring == 0        -0.039657   0.012794  -3.100  0.00493


incomp.aov<-lm(prop~relationship,data=incomp.dat)
anova(incomp.aov)
#Response: prop
#                Df Sum Sq   Mean Sq F value  Pr(>F)  
#relationship     2   0.04 0.0193527  3.7839 0.02274 *
#Residuals    73533 376.08 0.0051145  
incomp.tukey<-glht(incomp.aov, linfct = mcp(relationship="Tukey"))
summary(incomp.tukey)
#Linear Hypotheses:
#                                          Estimate Std. Error t value Pr(>|t|)
#mother-offspring - father-offspring == 0 -0.004177   0.015125  -0.276   0.9559
#unrelated - father-offspring == 0         0.014920   0.006278   2.377   0.0408
#unrelated - mother-offspring == 0         0.019097   0.013766   1.387   0.3270
                                 