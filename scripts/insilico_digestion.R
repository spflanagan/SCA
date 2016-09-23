#Author: Sarah P. Flanagan
#Last Updated: 23 September 2016
#Date Started: 23 September 2016
#Purpose: analyze in-silico digestion results

rm(list=ls())
setwd("~/Projects/SCA/results/insilico")

sd<-read.delim("SSC_sdigested.txt",header=T)
sd$RS<-paste(sd$Chrom,sd$PstIPos,sep=".")
dd<-read.delim("SSC_ddigested.txt",header=T)
dd$RS<-paste(dd$Chrom,dd$PstIPos,sep=".")
hist(sd$FragmentLength)

barplot(dpois(seq(0,10,1),lambda=0.000625))
