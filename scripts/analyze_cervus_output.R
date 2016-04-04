#Author: Sarah P. Flanagan
#Date: 15 February 2016
#Purpose: Analyze the output from CERVUS when it was run with various numbers
#of loci

rm(list=ls())
setwd("E:/ubuntushare/SCA/results/parentage")
maternity.files<-list.files(pattern="maternity.csv")
stats.files<-list.files(pattern="maternity.txt")

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

png("CervusStats.png",height=4,width=8,res=300, units="in")
par(mfrow=c(1,2),oma=c(2,2,2,2),mar=c(3,4,1,1))
boxplot(stats$Delta~stats$NumLoci,las=1,ylab="",xlab="")
mtext("Delta Cutoff",2,outer=F, line=2)
plot(stats$NumLoci, stats$AssignmentRate, xaxt='n',las=1,
	ylab="",xlab="", pch=19)
axis(1, at=c(50,100,200,400,800,1600))
mtext("Assignment Rate (%)",2,outer=F,line=2)
mtext("Number of Loci", 1,outer=T, line=-1)
dev.off()

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
