#Author: Sarah P. Flanagan
#Date: 16 Feb 2016
#Purpose: Analyze relatedness data - run band_sharing program first

setwd("B:/ubuntushare/SCA/results/relatedness")
rout<-read.table("genotypes99_10loci.rout.txt",sep="\t",header=T)

rout[rout$Ind1=="OFF083" & rout$Ind2=="PRM083",]


info<-read.table("genotypes99_10loci.txt",sep="\t",header=T)
dadkid<-read.table("../dad.kid.pairs.txt",sep="\t",header=F)

bandsharing<-function(genotypes, pairs){
	dki.list<-apply(pairs,1,function(x){
		dadinfo<-info[x[1] == genotypes$ID,]
		kidinfo<-info[x[2] == genotypes$ID,]
		return(rbind(dadinfo,kidinfo))
	})

	matches<-lapply(dki.list,function(y) {
		apply(y[,-1],2,function(x){ match<-x[1] == x[2] }) })
	matchcount<-unlist(lapply(matches, function(x) { 
		count<-length(x[x==TRUE]) }))
	shared<-as.numeric(matchcount)/20
	return(shared)
}

shared1600<-bandsharing(info,dadkid)

polyshared<-read.table("PolymorphicIn90PercInds.bandsharing.txt",header=T)
polyhweshared<-read.table("PolymorphicIn99PercIndsHWE.bandsharing.txt",header=T)

png("BandSharing.png",width=7,height=7,units="in",res=300)
par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
hist(polyshared$Shared, axes=F,main="",ylab="",xlab="",
	breaks=seq(-0.05,1.05,0.1),ylim=c(0,150))
axis(1,pos=0, xlim=c(-0.05,1.05))
axis(2,pos=-0.05,las=1,ylim=c(0,150))
mtext("Frequency",2,outer=T)

hist(polyshared$Incompatible, axes=F,main="",ylab="",xlab="",
	breaks=seq(-0.05,1.05,0.1),ylim=c(0,150))
axis(1,pos=0, xlim=c(-0.05,1.05))
axis(2,pos=-0.05,las=1,ylim=c(0,150))

hist(polyhweshared$Shared, axes=F,main="",ylab="",xlab="",
	breaks=seq(-0.05,1.05,0.1),ylim=c(0,150))
axis(1,pos=0, xlim=c(-0.05,1.05))
axis(2,pos=-0.05,las=1,ylim=c(0,150))
mtext("Band Sharing",1,outer=F,line=2)

hist(polyhweshared$Incompatible, axes=F,main="",ylab="",xlab="",
	breaks=seq(-0.05,1.05,0.1),ylim=c(0,150),xlim=c(-0.05,1))
axis(1,pos=0,xlim=c(-0.05,1.05))
axis(2,pos=-0.05,las=1,ylim=c(0,150))
mtext("Proportion Incompatible Loci", 1,outer=F,line=2)
dev.off()

