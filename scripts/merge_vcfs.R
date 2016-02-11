#Author: Sarah P. Flanagan
#Date: 9 February 2016
#Purpose: Merge FORMAT files exported from vcf.

setwd("E:/ubuntushare/SCA/results/biallelic")
vcf1<-read.delim("out.GT.FORMAT")
dups<-grep("align.1",colnames(vcf1))
vcf1<-vcf1[,-dups]
vcf2<-read.delim("fem.GT.FORMAT")
#vcf1$index<-paste(vcf1$CHROM,vcf1$POS,sep=".")
#vcf2$index<-paste(vcf2$CHROM,".vcf2$POS,sep=".")
vcf<-merge(vcf1,vcf2, by=c("CHROM","POS"))
vcf$index<-apply(vcf,1,function(x){ idx<-paste(x[1],".",x[2],sep="") })
addedon<-vcf[duplicated(vcf$index),"index"]
vcf<-vcf[!(vcf$index %in% addedon),]

write.table(vcf[,1:(ncol(vcf)-1)],"biallelic.vcf",col.names=T,row.names=F,
	quote=F,sep='\t')

