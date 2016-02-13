#Author: Sarah P. Flanagan
#Date: 12 February 2016
#Purpose: To generate an empirical allele frequency spectrum for all
#sequenced individuals from the output of gwsca_biallelic_vcf
rm(list=ls())

setwd("E:/ubuntushare/SCA/results/biallelic/")
gw.sum<-read.delim("gwsca_summary.txt")
sum.seq<-gw.sum[gw.sum$Pop == "ADULT" | gw.sum$Pop == "OFF",]
sum.seq<-sum.seq[!is.na(sum.seq$Hs),]
sum.seq$Count1<-sum.seq$Allele1Freq*sum.seq$N*2
overall.freq<-tapply(sum.seq$Count1,sum.seq$Locus,sum) 
overall.n<-tapply(sum.seq$N,sum.seq$Locus,sum) 
overall<-data.frame(Locus=names(overall.freq),Count1=overall.freq, 
	N=overall.n)
overall$Freq1<-overall$Count1/(2*overall$N)
overall$Freq2<-1-overall$Freq1		
max.af<-apply(overall[,c("Freq1","Freq2")],1,max)
af.hist<-hist(max.af, breaks=100)
af.export<-data.frame(breaks=af.hist$breaks[2:101], density=af.hist$density)
write.table(af.export, "empirical_allelefreqs.txt", 
	quote=F, col.names=F, row.names=F)