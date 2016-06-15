#Author: Sarah P. Flanagan
#Last Updated: 3 June 2016
#Date: 15 February 2016
#Purpose: Analyze the output from CERVUS when it was run with various numbers
#of loci

rm(list=ls())
setwd("E:/ubuntushare/SCA/results/parentage")

##############################################################################
#####CERVUS ANALYSIS
##############################################################################
stats.files<-list.files(pattern="\\d+_maternity.txt")

stats<-data.frame(NumLoci=numeric(),ConfidenceLevel=numeric(), Delta=numeric(),
	NumAssignments=numeric(), AssignmentRate=numeric())
for(i in 1: length(stats.files)){
	dat<- readLines(stats.files[i])
	num.loci<-gsub("gen(\\d+)_\\d+.maternity.txt","\\1",stats.files[i])
	start<-match("Mother given known father:", dat)
	#pull out info for strict only
	info<-unlist(strsplit(dat[(start+4)],"\\s+"))
	info<-unlist(strsplit(gsub("[()%]","",info),"[[:space:]]"))
	stats[i,]<-cbind(as.numeric(num.loci),as.numeric(info[2]),
		as.numeric(info[3]),as.numeric(info[4]), as.numeric(info[6]))
	rownames(stats)[i]<-stats.files[i]
}

#####PLOT CERVUS INFO
png("CervusStats.png",height=5,width=10,res=300, units="in")
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(3,3,1,0.2))
boxplot(stats$Delta~stats$NumLoci,las=1,ylab="",xlab="",xaxt='n')
axis(1,at=1:8,labels=c(50,100,150,200,300,400,800,1600),las=2)
#text(seq(0.8,7.8,1), par("usr")[1]-0.65, 
#	labels=c(50,100,150,200,300,400,800,1600), 
#	srt=35, pos=1, xpd=TRUE,tck=0.1)
mtext("LOD Cutoff",2,outer=F, line=2)
plot(stats$NumLoci, stats$AssignmentRate, xaxt='n',las=1,
	ylab="",xlab="", pch=19)
axis(1, at=c(50,100,150,200,300,400,800,1600),las=2)
lines(x=c(15,85),y=c(6.9,6.9),lwd=2)
lines(x=c(65,135),y=c(17.4,17.4),lwd=2)
lines(x=c(115,185),y=c(20.8,20.8),lwd=2)
lines(x=c(165,235),y=c(21.6,21.6),lwd=2)
lines(x=c(265,335),y=c(20.6,20.6),lwd=2)
lines(x=c(365,435),y=c(20.8,20.8),lwd=2)
lines(x=c(765,835),y=c(20.8,20.8),lwd=2)
lines(x=c(1565,1635),y=c(21,21),lwd=2)
mtext("Assignment Rate (%)",2,outer=F,line=2)
mtext("Number of Loci", 1,outer=T)
dev.off()

#####CHECK ALLELE FREQS
plot(stats$NumLoci, stats$AssignmentRate, xaxt='n',las=1,
	ylab="",xlab="",type='n')
stats$setID<-gsub("gen\\d+_(\\d+)_maternity.txt","\\1",rownames(stats))
text(labels=stats$setID,x=stats$NumLoci,y=stats$AssignmentRate)

af.files<-list.files(pattern="\\d+_allelefreqs.txt")
af.dat<-list()
for(i in 1:length(af.files)){
	nloci<-as.numeric(gsub("gen(\\d+)_\\d+_allelefreqs.txt","\\1",af.files[i]))
	af<-read.table(af.files[i],skip=10,header=T,nrows=nloci)
	af.dat[[i]]<-as.data.frame(af)
}
afsum<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ cbind(mean(x$k),mean(x$HObs)) })))
colnames(afsum)<-c("MeanK","MeanObsHet")
rownames(afsum)<-af.files
afsum$NumLoci<-as.numeric(
	gsub("gen(\\d+)_\\d+_allelefreqs.txt","\\1",rownames(afsum)))
afsum$setID<-gsub("gen\\d+_(\\d+)_allelefreqs.txt","\\1",rownames(afsum))
plot(afsum$NumLoci, afsum$MeanObsHet, xaxt='n',las=1,
	ylab="",xlab="",type='n',ylim=c(0,0.26))
text(labels=afsum$setID,x=afsum$NumLoci,y=afsum$MeanObsHet)
text(labels=stats$setID,x=stats$NumLoci,y=stats$AssignmentRate/100,col="blue")
axis(1, at=c(50,100,150,200,300,400,800,1600),las=2)

nalleles<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ summary(x$k) })))
rownames(nalleles)<-af.files
nalleles$NumLoci<-as.numeric(
	gsub("gen(\\d+)_\\d+_allelefreqs.txt","\\1",rownames(nalleles)))
nalleles$setID<-gsub("gen\\d+_(\\d+)_allelefreqs.txt","\\1",rownames(nalleles))
tapply(nalleles$Mean,nalleles$NumLoci,summary)

obshet<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ summary(x$HObs) })))
rownames(obshet)<-af.files
obshet$NumLoci<-as.numeric(
	gsub("gen(\\d+)_\\d+_allelefreqs.txt","\\1",rownames(obshet)))
obshet$setID<-gsub("gen\\d+_(\\d+)_allelefreqs.txt","\\1",rownames(obshet))
tapply(obshet$Mean,obshet$NumLoci,summary)

#####FULL DATASET
full.af<-read.table("PolymorhpicIn99PercIndsHWE_afreqs.txt",
	skip=10,header=T,nrows=1642)

full.dat<-read.table("PolymorphicIn99PercIndsHWE.txt",header=T)
ids<-as.character(full.dat$ID)
pairs<-data.frame(id1=character(),id2=character())
for(i in 1:(length(ids)-1)){
	ida<-rep(ids[i],(length(ids)-i))
	idb<-ids[(i+1):length(ids)]
	pairs<-rbind(pairs,cbind(ida,idb))
}
write.table(pairs,"pairwise.combinations.txt",col.names=F,row.names=F,
	quote=F,sep='\t')
##############################################################################
####OTHER ANALYSES
##############################################################################
maternity.files<-list.files(pattern="maternity.csv")
maternity.dat<-data.frame()
for(i in 1:length(maternity.files)){
	dat<-read.csv(maternity.files[i],skip=1,header=F)
	colnames(dat)<-scan(maternity.files[i],nlines=1,sep=",",
		what="character")
	sig<-dat[dat$"Trio confidence"=="*",]
	sig<-sig[,c("Offspring ID", "Candidate mother ID")]
	sig$NumLoci<-gsub("gen(\\d+)_\\d+.maternity.csv","\\1",maternity.files[i])
	maternity.dat<-rbind(maternity.dat, sig)
}
mat.split<-split(maternity.dat, factor(maternity.dat$"Candidate mother ID"))
summ.dat<- do.call("rbind", lapply(mat.split,function(x){ 
	sum<-summary(factor(x$"Offspring ID")) 
	min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ min(as.numeric(x$NumLoci)) }))
	max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ max(as.numeric(x$NumLoci)) }))
	y<-data.frame(MomID=levels(factor(x[,2])),OffID=names(sum),
		NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc)
	rownames(y)<-NULL
	return(y)
}))
rownames(summ.dat)<-NULL
write.csv(summ.dat, "Cervus_Summary.csv")

##split by locus instead of females
mat.loc<-split(maternity.dat, factor(maternity.dat$NumLoci))

mat.loc.sum<-lapply(mat.loc,function(ldf){ 
	df<-split(ldf,factor(ldf$"Candidate mother ID"))
	do.call("rbind", lapply(df, function(x){
		sum<-summary(factor(x$"Offspring ID")) 
		min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ min(as.numeric(x$NumLoci)) }))
		max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ max(as.numeric(x$NumLoci)) }))
		y<-data.frame(MomID=levels(factor(x[,2])),OffID=names(sum),
			NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc)
		rownames(y)<-NULL
		return(y)
	}))
})

png("CervusMaternitySummary.png",height=5,width=7,units="in",res=300)
par(mfrow=c(2,3), mar=c(2,2,2,2),oma=c(2,2,2,2))
hist(mat.loc.sum$`100`$NumAssignments, breaks=seq(0,11,0.5),xaxt='n',yaxt='n',
	main="",xlab="",ylab="",xlim=c(0,11),ylim=c(0,10))
legend("top",c("100 Loci"),bty='n')
mtext("Frequency", 2, outer=F,line=2,cex=0.75)
axis(1,pos=0)
axis(2, pos=0,las=1)
hist(mat.loc.sum$`200`$NumAssignments, breaks=seq(0,11,0.5),xaxt='n',yaxt='n',
	main="",xlab="",ylab="",xlim=c(0,11))
axis(1,pos=0)
axis(2, pos=0,las=1)
legend("top",c("200 Loci"),bty='n')
mtext("Number of Assignments", 1, outer=F,line=2,cex=0.75)
hist(mat.loc.sum$`1600`$NumAssignments, breaks=seq(0,11,0.5),xaxt='n',yaxt='n',
	main="",xlab="",ylab="",xlim=c(0,11))
axis(1,pos=0)
axis(2, pos=0,las=1)
legend("top",c("1600 Loci"),bty='n')

hist(summary(mat.loc.sum$`100`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30))
mtext("Frequency",2,outer=F,line=2,cex=0.75)
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
hist(summary(mat.loc.sum$`200`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30))
mtext("Inferred Mating Success",1,outer=F,line=2,cex=0.75)
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
hist(summary(mat.loc.sum$`1600`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30))
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
dev.off()




###Old graph, don't use
#png("CervusMaternitySummary.png",height=7,width=7,units="in",res=300)
#par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(3,3,1,1))
#boxplot(summ.dat$NumAssignments~summ.dat$MaxLoci,las=1,ylab="",xlab="")
#mtext("Number of Assignments",2,outer=F,line=2,cex=0.75)
#mtext("Maximum Number of Loci",1,outer=F, line=2,cex=0.75)

#boxplot(summ.dat$NumAssignments~summ.dat$MinLoci,las=1,ylab="",xlab="")
#mtext("Minimum Number of Loci",1,outer=F, line=2,cex=0.75)

#hist(summary(summ.dat$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),xaxt='n',yaxt='n',
#	main="",xlab="",ylab="")
#mtext("Number of Offspring Assigned",1,outer=F,line=2,cex=0.75)
#mtext("Frequency",2,outer=F,line=2,cex=0.75)
#axis(1,pos=0,at=c(1,2,3,4))
#axis(2,pos=0.5, las=1)

#hist(summ.dat$NumAssignments, breaks=20,xaxt='n',yaxt='n',
#	main="",xlab="",ylab="")
#mtext("Number of Assignments", 1, outer=F,line=2,cex=0.75)
#axis(1,pos=0)
#axis(2, pos=0,las=1)
#dev.off()
