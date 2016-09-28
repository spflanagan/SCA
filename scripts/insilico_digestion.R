#Author: Sarah P. Flanagan
#Last Updated: 28 September 2016
#Date Started: 23 September 2016
#Purpose: analyze in-silico digestion results

rm(list=ls())
setwd("~/Projects/SCA/results/insilico")
setwd("//VBOXSVR/Home/Projects/SCA/results/insilico")
source("../../scripts/plotting_functions.R")

sd<-read.delim("SSC_sdigested.txt",header=T)
sd$RS<-paste(sd$Chrom,sd$PstIPos,sep=".")
dd<-read.delim("SSC_ddigested.txt",header=T)
dd$RS<-paste(dd$Chrom,dd$PstIPos,sep=".")
hist(sd$FragmentLength)

barplot(dpois(seq(0,10,1),lambda=0.000625))

#analyze the fsts

summ<-data.frame(read.table("ssc_insilico_summstats.txt",sep='\t',
	header=T,stringsAsFactors=F))
summ$Fst<-as.numeric(summ$Fst)
ifsts<-summ[!is.na(summ$Fst),]
ifsts<-ifsts[as.numeric(ifsts$Fst)>=0,]
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
scaffs<-levels(as.factor(ifsts[,1]))
scaffs[1:22]<-lgs

isf<-fst.plot(ifsts,ci.dat=c(-10,10),sig.col=c("black","black"),fst.name="Fst",
	chrom.name="Chrom",bp.name="Pos",axis.size=1,
	groups=as.factor(scaffs[scaffs %in% levels(factor(ifsts[,1]))]))
