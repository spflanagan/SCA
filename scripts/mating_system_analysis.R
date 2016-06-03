#Author: Sarah P. Flanagan
#Date: 1 June 2016
#Last updated: 1 June 2016
#Purpose: Analyze mating system data from single population of scovelli

rm(list=ls())
setwd("../results/parentage")
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
par(mfrow=c(2,2),mar=c(2,2,2,0),oma=c(2,2,0.5,0.5))
#plot males
plot(mal.dat$NumMates,mal.dat$No.Offspring,pch=17,xlim=c(0,3),ylim=c(-1,70),
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.25,y=65,"Males")
clip(0,1,0,70)
abline(lm(No.Offspring~NumMates,dat=mal.dat))
mtext("Number of Offspring",2,outer=F,line=1.5,cex=0.85)
mtext("Number of Mates",1,outer=F,line=1,cex=0.85)

#plot females with no zeros
plot(fem.dat$NumMates,fem.dat$No.Offspring,xlim=c(0,3),ylim=c(-1,70),pch=19,
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.25,y=65,"Females")
clip(1,3,0,70)
#abline(a=0,b=17.69) #use the Bateman gradient slope? 
abline(lm(No.Offspring~NumMates,dat=fem.dat))
clip(0,3,0,70)
abline(lm(No.Offspring~NumMates,dat=femz.dat),lty=2,col="blue")

legend(x=1.9,y=10,lty=c(1,2),col=c("black","blue"),
	c("Without Zeroes","With Zeroes"),bty='n',ncol=1)

mtext("Number of Mates",1,outer=F,line=1,cex=0.85)

#dev.off()
#########PLOT HISTOGRAMS
#png("MS_RS_distributions.png",height=10,width=5,units="in",res=300)
#par(mfrow=c(2,1),mar=c(2,2,2,2),oma=c(2,2,0.5,0.5))
#plot observed female mating success
hist(femz.dat$NumMates,breaks=seq(-0.5,3.5,1),axes=F,main="",xlab="",ylab="",
	ylim=c(0,41))
axis(1,at=seq(0,3,1),lab=seq(0,3,1),pos=0)
axis(2,pos=-0.5,las=1,at=seq(0,40,10),lab=seq(0,40,10))
mtext("Number of Mates",1,outer=F,line=1.5,cex=0.85)
mtext("Number of Females",2,outer=F,line=1.5,cex=0.85)

#plot observed female reproductive mating success
hist(femz.dat$No.Offspring,axes=F,main="",xlab="",ylab="",ylim=c(0,41))
axis(1,pos=0)
axis(2,pos=0,las=1,at=seq(0,40,10),lab=seq(0,40,10))
mtext("Number of Offspring",1.5,outer=F,line=2,cex=0.85)

dev.off()



#############################################################################
#########BAND SHARING
#############################################################################
setwd("../relatedness/")
all.bs<-read.delim("PolymorphicIn99PercIndsHWE.allcombos.bandsharing.txt")
fo.bs<-read.delim("PolymorphicIn99PercIndsHWE.bandsharing.txt")
all.bs$combo<-paste(all.bs$Father,all.bs$Offspring,sep=".")
all.bs$combo2<-paste(all.bs$Offspring,all.bs$Father,sep=".")

fo.bs$combo<-paste(fo.bs$Father,fo.bs$Offspring,sep=".")

unr.bs<-all.bs[!(all.bs$combo %in% fo.bs$combo) & 
	!(all.bs$combo2 %in% fo.bs$combo),]
unr.bs$combo<-rep("unrelated",nrow(unr.bs))
fo.bs$combo<-rep("father-offspring",nrow(fo.bs))
all.bs<-rbind(unr.bs[,colnames(unr.bs)[colnames(unr.bs)!="combo2"]],fo.bs)
all.bs[all.bs$combo=="father-offspring" & all.bs$Shared < 0.2,
	"combo"]<-"unrelated"

boxplot(all.bs$Shared~all.bs$combo)
boxplot(all.bs$Incompatible~all.bs$combo)

shared.dat<-data.frame(prop=all.bs$Shared,relationship=all.bs$combo,
	type="shared")
incomp.dat<-data.frame(prop=all.bs$Incompatible,relationship=all.bs$combo,
	type="incompatible")
all.dat<-rbind(shared.dat,incomp.dat)
boxplot(all.dat$prop)
require(ggplot2)
ggplot(data = all.dat, aes(x=relationship, y=prop)) + 
	geom_boxplot(aes(fill=type))
