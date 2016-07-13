#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Date: 1 July 2016
#Purpose: QC of RAD-seq data
#check heterozygosity vs read depth
#distribution of reads per allele in heterozygotes
#calculate coverage per locus

setwd("E:/ubuntushare/SCA/results/quality_control/")
##FILES
hwe.map<-read.delim("../stacks/plink.map",header=F)
colnames(hwe.map)<-c("Chrom","Locus","Dist","BP")
#hwe.map$Locus<-gsub("(\\d+)_(\\d+)","\\1.\\2",hwe.map$Locus)
hwe.map$LocID<-paste("LG",hwe.map$Chrom,".",hwe.map$Locus,sep="")
sumstats<-read.delim("../stacks/batch_1.sumstats.tsv",comment.char="#",header=F)
colnames(sumstats)<-c("BatchID","Locus","Chr","BP","Col","PopID","PNuc",
	"QNuc","N","P","ObsHet","ObsHom","ExpHet","ExpHom","Pi","SmoothedPi",
	"SmoothedPiP-value","Fis","SmoothedFis","SmoothedFisP-value","Private")
sumstats$PlinkID<-paste(sumstats$Locus,sumstats$Col,sep="_")
sumstats$LocusID<-paste(sumstats$Chr,sumstats$Locus,sumstats$BP,sep=".")
vcf<-read.delim("stacks/batch_1.vcf",header=F,comment.char="#",sep='\t')
header.start<-grep("#CHROM",scan("stacks/batch_1.vcf",what="character"))
header<-scan("stacks/batch_1.vcf",what="character")[header.start:
	(header.start+ncol(vcf)-1)]
colnames(vcf)<-header

lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
cov<-read.delim("coverage_info.txt",header=T,stringsAsFactors=F)
cov$LocID<-paste(cov$Chrom,cov$Locus,cov$BP,sep=".")



#hwe
hwe<-merge(sumstats[,c("PlinkID","LocusID")],hwe.map,by.x="PlinkID",
	by.y="Locus")
cov.hwe<-cov[cov$LocID %in% hwe$LocusID,]

#prune to keep only those found >100 adults and > 100 off
sum.prune<-sumstats[sumstats$PopID=="PRM" | sumstats$PopID=="OFF",]
sum.prune<-sum.prune[sum.prune$N > 100,]
sum.prune<-sum.prune[duplicated(sum.prune$LocusID),]

#Now prune based on allele frequencies
sum.afs<-tapply(sumstats$P,sumstats$LocusID,mean)
sum.afs<-sum.afs[sum.afs>0.05 & sum.afs<0.95]
sum.prune<-sum.prune[sum.prune$LocusID %in% names(sum.afs),]

#prune for coverage
cov.prune<-cov[cov$LocID %in% sum.prune$LocusID,]
cov.prune<-cov.prune[cov.prune$AvgDepthPerInd >=1 & 
	cov.prune$AvgDepthPerInd <=100,]

#get genotypic likelihoods
gl<-apply(vcf,c(1,2),gsub,pattern="./.:\\d+:.,.:.,.,.",replacement="NA")
gl1<-apply(gl,c(1,2),gsub,
	pattern="\\d+/\\d+:\\d+:\\d+,\\d+:(*)",
	replacement="\\1")
gl2<-data.frame(apply(gl1,c(1,2),gsub,
	pattern="\\.,(\\d+\\.\\d+),\\.",
	replacement="\\1"))
gl3<-apply(gl2,2,as.character)
gl3<-apply(gl3[,10:ncol(gl3)],2,as.numeric)#I'll remove those with LR < 10
means<-cbind(vcf[,1:3],rowMeans(gl3,na.rm=T))
keep.lr<-means[means[,4]>20,]
keep.lr$Locus<-paste(keep.lr$`#CHROM`,keep.lr$ID,keep.lr$POS,sep=".")
cov.prune<-cov.prune[cov.prune$LocID %in% keep.lr$Locus,]#22682


##gentoypes
gt<-apply(gl,c(1,2),gsub,
	pattern="(\\d+/\\d+):\\d+:\\d+,\\d+:(.*)",
	replacement="\\1")
ad<-apply(gl,c(1,2),gsub,
	pattern="\\d+/\\d+:\\d+:(\\d+,\\d+):(.*)",
	replacement="\\1")

gt.sub<-data.frame(gt[gt[,3] %in% sample(gt[,3], 100,replace=F),])
ad.sub<-data.frame(ad[ad[,3] %in% gt.sub$ID,])

alleles<-list()
for(j in 1:nrow(gt.sub)){
	cab<-data.frame(allele0=numeric(),allele1=numeric(),ratio=numeric(),
		stringsAsFactors=F)	
	this.row<-data.frame(t(gt.sub[j,10:ncol(gt.sub)]))
	hets<-this.row[this.row[,1] %in% c("0/1","1/0"),]
	for(i in 1:length(hets)){
		index<-names(hets)[i]
		a<-substr(gt.sub[j,index],1,1)
		b<-substr(gt.sub[j,index],3,3)
		ca<-as.numeric(strsplit(as.character(ad.sub[j,index]),",")[[1]][1])	
		cb<-as.numeric(strsplit(as.character(ad.sub[j,index]),",")[[1]][2])	
		if(a == "0") {
			cab[i,1]<-ca
		} else {
			cab[i,2]<-ca
		}
		if(b == "0") {
			cab[i,1]<-cb
		} else {
			cab[i,2]<-cb
		}
		cab[i,"ratio"]<-ca/(ca+cb)
	} 
	alleles[[j]]<-as.vector(cab[,"ratio"])
}

all.in<-c(seq(1,length(alleles),49),length(alleles)+1)
for(j in 1:(length(all.in)-1)){
png(paste("HetRatios.",j,".png",sep=""),height=10.5,width=8,units="in",res=300)
par(mfrow=c(7,7),mar=c(1,1,1,1),oma=c(2,2,2,2))
for(i in all.in[j]:(all.in[j+1]-1)){
	h1<-hist(alleles[[i]],breaks=10,plot=F)
	plot(h1,col=alpha("slategray",0.5),border=alpha("slategray",0.5),
		ylim=c(0,max(h1$counts)),main="",xlab="",ylab="",
		axes=F,xlim=c(0,1))
	axis(1,pos=0,at=seq(0,1,0.5),
		labels=seq(0,1,0.5),xlim=c(0,1))
	axis(2,pos=0,las=1)
}
mtext("Number of Alleles",1,outer=T,line=1)
mtext("Number of Individuals",2,outer=T,line=1)
dev.off()
}


##Look at coverage vs. heterozygosity
plot(log(cov.prune$AvgDepthPerInd),cov.prune$PropHet,
	col=alpha(as.numeric(as.factor(cov.prune$Chrom)),0.5),pch=19)
#need to break it up
breaks<-c(seq(0,19,1),seq(20,100,10),seq(150,3000,50))
cov$DepthGroup<-cov$AvgDepthPerInd
depth.plot<-NULL
for(i in 1:(length(breaks)-1)){
	depth.plot<-rbind(depth.plot,
		c(mean(breaks[i:(i+1)]),
		mean(cov$PropHet[cov$DepthGroup > breaks[i] & 
			cov$DepthGroup < breaks[i+1]])))
}
png("Het_cov.png",height=7,width=7,units="in",res=300)
plot(depth.plot[,1],depth.plot[,2],
	xlab="Average Allele Depth Per Individual",
	ylab="Proportion Called Heterozygous",pch=19,las=1,lwd=1.3,cex=1.3)
dev.off()

##Check coverage per chromosome
cov.lg<-cov[cov$Chrom %in% lgs,]
png("Coverage_histograms.png",height=10,width=10,units="in",res=300)
par(mfrow=c(4,6),mar=c(2,2,2,2),oma=c(2,2,2,2),las=1,lwd=1.3)
tapply(log(cov.lg$TotalDepth), cov.lg$Chrom,hist, main="",xlab="",ylab="")
hist(log(cov.lg$TotalDepth),main="",xlab="",ylab="",col="dodgerblue")
hist(log(cov$TotalDepth),main="",xlab="",ylab="",col="slategrey")
mtext("log(Total Sequencing Depth)",1,outer=T,cex=0.75)
mtext("Number of Loci",2,outer=T,las=0,line=1.1,cex=0.75)
dev.off()

##Distribution of reads per allele in heterozygotes
png("ReadDepthRatioHets.png",height=7,width=7,units="in",res=300)
par(mar=c(2,2,2,2),oma=c(1,3,1,0),lwd=1.3,cex=1.3)
hist(log(cov$HetDepthRatio), main="",xlab="",ylab="",las=1,xaxt='n')
axis(1,pos=0)
mtext("Read Depth Ratio in Heterozygotes (log)",1,outer=F,line=1.5,cex=1.3)
mtext("Number of Heterozygous Loci",2,outer=F,line=3.5,cex=1.3)
legend("topright",ncol=1,bty='n',
	c("Mean = 1.288","Median=1.13", "Min=0","Max=184.71",">5=0.43%"))
dev.off()