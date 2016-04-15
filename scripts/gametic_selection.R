#Author: Sarah P. Flanagan
#Date: 15 April 2016
#Purpose: Calculate gametic selection

setwd("E:/ubuntushare/SCA/results/biallelic")
vcf<-read.table("biallelic_merge.vcf",header=T)
pairs.names<-read.table("../dad.kid.pairs.fullnames.txt")

gametic.selection<-function(vcf.file,parent.off.pair){
	print(parent.off.pair)
	these<-c(as.character(parent.off.pair[1]),as.character(parent.off.pair[2]))
	pair<-cbind(vcf.file[,1:9],vcf.file[,these])
	dad.het<-pair[pair[,10]=="0/1" | pair[,10]=="1/0",]
	dad.het.off.het<-dad.het[dad.het[,11]=="0/1" | dad.het[,11]=="1/0",]
	prop<-nrow(dad.het.off.het)/nrow(dad.het)
}

gs<-apply(pairs.names,1,gametic.selection,vcf.file=vcf)