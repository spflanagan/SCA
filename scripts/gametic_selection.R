#Author: Sarah P. Flanagan
#Date: 15 April 2016
#Purpose: Calculate gametic selection

setwd("E:/ubuntushare/SCA/results/biallelic")
vcf<-read.table("biallelic_merge.vcf",header=T)
pairs.names<-read.table("../dad.kid.pairs.fullnames.txt",stringsAsFactors=F)
pairs.names$V1[pairs.names$V1=="sample_PRM086-23_align"]<-"sample_PRM086.23_align"
remove<-"sample_PRM077_align"
pairs.names<-pairs.names[pairs.names$V1 != remove,]
prg.n.loc<-read.table("LociInPregMales.txt",header=T)
juv.n.loc<-read.table("LociInOffspring.txt",header=T)
prg.n.loc$Locus<-paste(prg.n.loc$Chrom,prg.n.loc$LocID,prg.n.loc$Pos,sep=".")
juv.n.loc$Locus<-paste(juv.n.loc$Chrom,juv.n.loc$LocID,juv.n.loc$Pos,sep=".")
vcf$Locus<-paste(vcf$CHROM,vcf$ID,vcf$POS,sep=".")

vcf.prune<-vcf[vcf$Locus %in% prg.n.loc$Locus & 
	vcf$Locus %in% juv.n.loc$Locus,]
vcf.prune<-vcf.prune[,colnames(vcf.prune)!="Locus"]

gametic.selection<-function(vcf.file,parent.off.pair){
	print(parent.off.pair)
	these<-c(as.character(parent.off.pair[1]),as.character(parent.off.pair[2]))
	pair<-cbind(vcf.file[,1:9],vcf.file[,these])
	dad.het<-pair[pair[,10]=="0/1" | pair[,10]=="1/0",]
	dad.het.off.het<-dad.het[dad.het[,11]=="0/1" | dad.het[,11]=="1/0",]
	prop<-nrow(dad.het.off.het)/nrow(dad.het)
}

gs<-apply(pairs.names,1,gametic.selection,vcf.file=vcf.prune)
t.test(gs,mu=0.5)
sem<-sd(gs)/sqrt(length(gs))