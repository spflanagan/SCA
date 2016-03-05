#Author: Sarah P. Flanagan
#Date: 29 February 2016
#Purpose: Troubleshooting weird patterns in gwsca_biallelic_analysis

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")

#compare the two datasets
gw.ds.fst<-read.delim("gwsca_fsts_datasets.txt")
gw.ds.sum<-read.delim("gwsca_summary_datasets.txt")
gw.ds.fst$Locus<-paste(gw.ds.fst$Chrom,gw.ds.fst$Pos,sep=".")
gw.ds.sum$Locus<-paste(gw.ds.sum$Chrom,gw.ds.sum$Pos,sep=".")
gw.ds.sum<-gw.ds.sum[!is.na(gw.ds.sum$AA),]
gw.ds.sum$AAexp<-gw.ds.sum$Allele1Freq*gw.ds.sum$Allele1Freq
gw.ds.sum$aaexp<-gw.ds.sum$Allele2Freq*gw.ds.sum$Allele2Freq
gw.ds.sum$Aaexp<-1-gw.ds.sum$aaexp-gw.ds.sum$AAexp
gw.ds.sum$chi<-(((gw.ds.sum$AA-gw.ds.sum$AAexp)^2)/gw.ds.sum$AAexp)+
	(((gw.ds.sum$Aa-gw.ds.sum$Aaexp)^2)/gw.ds.sum$Aaexp)+
	(((gw.ds.sum$aa-gw.ds.sum$aaexp)^2)/gw.ds.sum$aaexp)
gw.ds.sum$chi.result<-1-pchisq(gw.ds.sum$chi,1) #biallelic, df=1
sum.sum<-tapply(gw.ds.sum$N,gw.ds.sum$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 87]#total N is 58*2<-this is 75%
psum.prune<-gw.ds.sum[gw.ds.sum$Locus %in% names(sum.sum),]
psum.list<-split(psum.prune, psum.prune$Pop)
dr.n<-psum.list$ddRAD[psum.list$ddRAD$N > 383,]
or.n<-psum.list$oRAD[psum.list$oRAD$N>58,]
#Now prune based on allele frequencies
dr.n<-dr.n[dr.n$Allele1Freq > 0.05 & dr.n$Allele1Freq < 0.95,]
or.n<-or.n[or.n$Allele1Freq > 0.05 & or.n$Allele1Freq < 0.95,]
#comparisons
ds.prune<-gw.ds.fst[gw.ds.fst$Locus %in% dr.n$Locus &
	gw.ds.fst$Locus %in% or.n$Locus, ]
ds.prune<-ds.prune[ds.prune$ddRAD.oRAD>0,]

outliers<-ds.prune[ds.prune$ddRAD.oRAD > 0.3,"Locus"]
out.sum<-gw.ds.sum[gw.ds.sum$Locus %in% outliers,]
out.sum<-out.sum[out.sum$Pop == "ddRAD" | out.sum$Pop == "oRAD",]
ho.diff<-tapply(out.sum$Ho,out.sum$Locus,function(x){ diff<-x[2]-x[1] })
hs.diff<-tapply(out.sum$Hs,out.sum$Locus,function(x){ diff<-x[2]-x[1] })

png("ddRADvoRAD.png",height=7,width=7,units="in",res=300)
plot(ds.prune$ddRAD.oRAD,xlab="",ylab="Fst")
dev.off()

#pstI only
gw.psti.fst<-read.delim("gwsca_fsts_psti.txt")
gw.psti.sum<-read.delim("gwsca_summary_psti.txt")
gw.psti.fst$Locus<-paste(gw.psti.fst$Chrom,gw.psti.fst$Pos,sep=".")
gw.psti.sum$Locus<-paste(gw.psti.sum$Chrom,gw.psti.sum$Pos,sep=".")
gw.psti.sum<-gw.psti.sum[!is.na(gw.psti.sum$AA),]
gw.psti.sum$AAexp<-gw.psti.sum$Allele1Freq*gw.psti.sum$Allele1Freq
gw.psti.sum$aaexp<-gw.psti.sum$Allele2Freq*gw.psti.sum$Allele2Freq
gw.psti.sum$Aaexp<-1-gw.psti.sum$aaexp-gw.psti.sum$AAexp
gw.psti.sum$chi<-(((gw.psti.sum$AA-gw.psti.sum$AAexp)^2)/gw.psti.sum$AAexp)+
	(((gw.psti.sum$Aa-gw.psti.sum$Aaexp)^2)/gw.psti.sum$Aaexp)+
	(((gw.psti.sum$aa-gw.psti.sum$aaexp)^2)/gw.psti.sum$aaexp)
gw.psti.sum$chi.result<-1-pchisq(gw.psti.sum$chi,1) #biallelic, df=1
sum.sum<-tapply(gw.psti.sum$N,gw.psti.sum$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 87]#total N is 58*2<-this is 75%
psum.prune<-gw.psti.sum[gw.psti.sum$Locus %in% names(sum.sum),]
psum.list<-split(psum.prune, psum.prune$Pop)
pfem.n<-psum.list$FEM[psum.list$FEM$N > 30,]
pmal.n<-psum.list$PRM[psum.list$PRM$N>28,]
#Now prune based on allele frequencies
pfem.n<-pfem.n[pfem.n$Allele1Freq > 0.05 & pfem.n$Allele1Freq < 0.95,]
pmal.n<-pmal.n[pmal.n$Allele1Freq > 0.05 & pmal.n$Allele1Freq < 0.95,]
#comparisons
pfm.prune<-gw.psti.fst[gw.psti.fst$Locus %in% pmal.n$Locus &
	gw.psti.fst$Locus %in% pfem.n$Locus, ]
pfm.prune<-pfm.prune[pfm.prune$FEM.PRM>0,]


#all ddRAD inds
dgw.fst<-read.delim("gwsca_fsts_ddRAD.txt")
dgw.sum<-read.delim("gwsca_summary_ddRAD.txt")
dgw.fst$Locus<-paste(dgw.fst$Chrom,dgw.fst$Pos,sep=".")
dgw.sum$Locus<-paste(dgw.sum$Chrom,dgw.sum$Pos,sep=".")
dgw.sum<-dgw.sum[!is.na(dgw.sum$AA),]
dgw.sum$AAexp<-dgw.sum$Allele1Freq*dgw.sum$Allele1Freq
dgw.sum$aaexp<-dgw.sum$Allele2Freq*dgw.sum$Allele2Freq
dgw.sum$Aaexp<-1-dgw.sum$aaexp-dgw.sum$AAexp
dgw.sum$chi<-(((dgw.sum$AA-dgw.sum$AAexp)^2)/dgw.sum$AAexp)+
	(((dgw.sum$Aa-dgw.sum$Aaexp)^2)/dgw.sum$Aaexp)+
	(((dgw.sum$aa-dgw.sum$aaexp)^2)/dgw.sum$aaexp)
dgw.sum$chi.result<-1-pchisq(dgw.sum$chi,1) #biallelic, df=1
sum.sum<-tapply(dgw.sum$N,dgw.sum$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 87]#total N is 58*2<-this is 75%
dsum.prune<-dgw.sum[dgw.sum$Locus %in% names(sum.sum),]
dsum.list<-split(dsum.prune, dsum.prune$Pop)
dfem.n<-dsum.list$FEM[dsum.list$FEM$N > 30,]
dmal.n<-dsum.list$PRM[dsum.list$PRM$N>28,]
#Now prune based on allele frequencies
dfem.n<-dfem.n[dfem.n$Allele1Freq > 0.05 & dfem.n$Allele1Freq < 0.95,]
dmal.n<-dmal.n[dmal.n$Allele1Freq > 0.05 & dmal.n$Allele1Freq < 0.95,]
#comparisons
dfm.prune<-dgw.fst[dgw.fst$Locus %in% dmal.n$Locus &
	dgw.fst$Locus %in% dfem.n$Locus, ]
dfm.prune<-dfm.prune[dfm.prune$FEM.PRM>0,]

#ddRAD subset
dsgw.fst<-read.delim("gwsca_fsts_ddRADsub.txt")
dsgw.sum<-read.delim("gwsca_summary_ddRADsub.txt")
dsgw.fst$Locus<-paste(dsgw.fst$Chrom,dsgw.fst$Pos,sep=".")
dsgw.sum$Locus<-paste(dsgw.sum$Chrom,dsgw.sum$Pos,sep=".")
dsgw.sum<-dsgw.sum[!is.na(dsgw.sum$AA),]
dsgw.sum$AAexp<-dsgw.sum$Allele1Freq*dsgw.sum$Allele1Freq
dsgw.sum$aaexp<-dsgw.sum$Allele2Freq*dsgw.sum$Allele2Freq
dsgw.sum$Aaexp<-1-dsgw.sum$aaexp-dsgw.sum$AAexp
dsgw.sum$chi<-(((dsgw.sum$AA-dsgw.sum$AAexp)^2)/dsgw.sum$AAexp)+
	(((dsgw.sum$Aa-dsgw.sum$Aaexp)^2)/dsgw.sum$Aaexp)+
	(((dsgw.sum$aa-dsgw.sum$aaexp)^2)/dsgw.sum$aaexp)
dsgw.sum$chi.result<-1-pchisq(dsgw.sum$chi,1) #biallelic, df=1
sum.sum<-tapply(dsgw.sum$N,dsgw.sum$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 87]#total N is 58*2<-this is 75%
dssum.prune<-dsgw.sum[dsgw.sum$Locus %in% names(sum.sum),]
dssum.list<-split(dssum.prune, dssum.prune$Pop)
dsfem.n<-dssum.list$FEM[dssum.list$FEM$N > 30,]
dsmal.n<-dssum.list$PRM[dssum.list$PRM$N>28,]
#Now prune based on allele frequencies
dsfem.n<-dsfem.n[dsfem.n$Allele1Freq > 0.05 & dsfem.n$Allele1Freq < 0.95,]
dsmal.n<-dsmal.n[dsmal.n$Allele1Freq > 0.05 & dsmal.n$Allele1Freq < 0.95,]
#comparisons
dsfm.prune<-dsgw.fst[dsgw.fst$Locus %in% dsmal.n$Locus &
	dsgw.fst$Locus %in% dsfem.n$Locus, ]
dsfm.prune<-dsfm.prune[dsfm.prune$FEM.PRM>0,]

png("ddRADvoRAD.png",height=7,width=21,units="in",res=300)
par(mfrow=c(1,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(pfm.prune$FEM.PRM,ylab="",xlab="",axes=F,ylim=c(0,0.6))
axis(1,pos=0)
axis(2,pos=0,las=2,ylim=c(0,0.6))
mtext("original RAD",3,outer=F)
abline(h=0.025,col="grey")
abline(h=0.05,col="grey")
plot(dfm.prune$FEM.PRM,ylab="",xlab="",axes=F,ylim=c(0,0.6))
axis(1,pos=0)
axis(2,pos=0,las=2,ylim=c(0,0.6))
mtext("ddRAD",3,outer=F)
abline(h=0.025,col="grey")
abline(h=0.05,col="grey")
plot(dsfm.prune$FEM.PRM,ylab="",xlab="",axes=F,ylim=c(0,0.6))
axis(1,pos=0)
axis(2,pos=0,las=2,ylim=c(0,0.6))
abline(h=0.025,col="grey")
abline(h=0.05,col="grey")
mtext("ddRAD Subset",3,outer=F)
mtext("Index",1,outer=T,line=2)
mtext("Fst",2,outer=T,line=2)
dev.off()


