#Author: Sarah P. Flanagan
#Date: 16 Feb 2016
#Purpose: Analyze relatedness data

setwd("E:/ubuntushare/SCA/results/relatedness")
rout<-read.table("genotypes99_10loci.rout.txt",sep="\t",header=T)

rout[rout$Ind1=="OFF083" & rout$Ind2=="PRM083",]
