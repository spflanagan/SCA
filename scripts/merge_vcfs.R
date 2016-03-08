#Author: Sarah P. Flanagan
#Date: 9 February 2016
#Purpose: Merge FORMAT files exported from vcf.

setwd("E:/ubuntushare/SCA/results/biallelic")
vcf1<-read.delim("biallelic_maternal.GT.FORMAT")
vcf2<-read.delim("fem.GT.FORMAT")
vcf1$index<-paste(vcf1$CHROM,vcf1$POS,sep=".")
vcf2$index<-paste(vcf2$CHROM,vcf2$POS,sep=".")
vcf<-merge(vcf1,vcf2, by="index")

#remove the oddly duplicated ones
addedon<-vcf[duplicated(vcf$index),"index"]
if(!is.null(dim(addedon))) vcf<-vcf[!(vcf$index %in% addedon),]

#otherwise, clean it up and write to file.
drops<-c("index","CHROM.y","POS.y")
vcf<-vcf[,!(colnames(vcf) %in% drops)]
colnames(vcf)[1:2]<-c("CHROM","POS")
write.table(vcf,"biallelic.gt.vcf",col.names=T,row.names=F,
	quote=F,sep='\t')


