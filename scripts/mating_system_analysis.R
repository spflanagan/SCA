#Author: Sarah P. Flanagan
#Date: 1 June 2016
#Last updated: 1 June 2016
#Purpose: Analyze mating system data from single population of scovelli

rm(list=ls())
setwd("B:/ubuntushare/SCA/results/parentage")
dat<-read.delim("batemanator_input.rerun.txt")
fem.dat<-dat[substr(dat$Fish.ID,1,3)=="FEM",]
mal.dat<-dat[substr(dat$Fish.ID,1,3)!="FEM",]
zdat<-read.delim("batemanator_input.rerun_with0s.txt")
femz.dat<-zdat[substr(zdat$Fish.ID,1,3)=="FEM",]
sim.dat<-read.delim("batemanater_output_est.txt")

#############################################################################
#########PLOT BATEMAN GRADIENTS
#############################################################################

par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
#plot males
plot(mal.dat$NumMates,mal.dat$No.Offspring,pch=17,xlim=c(0,3),ylim=c(-1,70),
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.1,y=65,"A")
mtext("Males",3,cex=0.75)
clip(0,1,0,70)
abline(lm(No.Offspring~NumMates,dat=mal.dat))

#plot females with no zeros
plot(fem.dat$NumMates,fem.dat$No.Offspring,xlim=c(0,3),ylim=c(-1,70),pch=19,
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.1,y=65,"B")
mtext("Females, raw data, with zeros",3,cex=0.75)
clip(1,3,0,70)
#abline(a=0,b=17.69) #use the Bateman gradient slope? 
abline(lm(No.Offspring~NumMates,dat=fem.dat))

#plot females with zeros
plot(femz.dat$NumMates,femz.dat$No.Offspring,xlim=c(0,3),ylim=c(-1,70),pch=19,
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.1,y=65,"C")
mtext("Females, raw data, with zeros",3,cex=0.75)
clip(0,3,0,70)
#abline(a=0,b=17.69) #use the Bateman gradient slope? 
abline(lm(No.Offspring~NumMates,dat=fem.dat))

#plot estimated data
plot(sim.dat$NumMates,sim.dat$No..Offspring,xlim=c(0,3),ylim=c(-1,70),pch=19,
	xlab="",ylab="",xaxt='n',yaxt='n',axes=F)
axis(1,at=c(0,1,2,3),c(0,1,2,3),las=1,pos=0)
axis(2,las=1,pos=0)
text(x=0.1,y=65,"D")
mtext("Females, estimated",3,cex=0.75)
clip(1,3,0,70)
#abline(a=0,b=17.69) #use the Bateman gradient slope? 
abline(lm(No.Offspring~NumMates,dat=fem.dat))



mtext("Number of Mates",1,outer=T,cex=0.75)
mtext("Number of Offspring",2,outer=T,cex=0.75)

#############################################################################
#########PLOT HISTOGRAMS
#############################################################################

par(mfrow=c(2,1),mar=c(2,2,2,2),oma=c(2,2,2,2))
#plot observed female mating success
hist(femz.dat$NumMates,breaks=seq(-0.5,3.5,1),axes=F,main="",xlab="",ylab="",
	ylim=c(0,41))
axis(1,at=seq(0,3,1),lab=seq(0,3,1),pos=0)
axis(2,pos=-0.5,las=1,at=seq(0,40,10),lab=seq(0,40,10))
mtext("Number of Mates",1,outer=F,line=2)

#plot observed female reproductive mating success
hist(femz.dat$No.Offspring,axes=F,main="",xlab="",ylab="",ylim=c(0,41))
axis(1,pos=0)
axis(2,pos=0,las=1,at=seq(0,40,10),lab=seq(0,40,10))
mtext("Number of Offspring",1,outer=F,line=2)

mtext("Number of Individuals",2,outer=T)



