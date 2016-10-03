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

summstats.plot<-function(file.list,make.png=T,file.name="")
{
  lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
         "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
         "LG22")
  if(make.png==T) png(file.name,height=10,width=7.5,units="in",res=300)
  par(mfrow=c(length(file.list),1),oma=c(2,2,2,2),mar=c(2,2,2,2))
  for(i in 1:length(file.list)){
    summ<-data.frame(read.table(file.list[i],sep='\t',
                                header=T,stringsAsFactors=F))
    summ$Fst<-as.numeric(summ$Fst)
    ifsts<-summ[!is.na(summ$Fst),]
    ifsts<-ifsts[as.numeric(ifsts$Fst)>=0,]
    
    scaffs<-levels(as.factor(ifsts[,1]))
    scaffs[1:22]<-lgs
    
    i1<-fst.plot(ifsts,ci.dat=c(-10,10),sig.col=c("black","black"),fst.name="Fst",
                 chrom.name="Chrom",bp.name="Pos",axis.size=1,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(ifsts[,1]))]))
    last<-0
    for(ii in 1:length(lgs)){
      text(x=mean(i1[i1$Chrom ==lgs[ii],"Pos"]),y=-0.01,
           labels=lgs[ii], adj=1, xpd=TRUE,srt=90,cex=1)
      last<-max(i1[i1$Chrom ==lgs[ii],"Pos"])
    }
    mtext(file.list[i],3)
  }
  mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
  if(make.png==T) dev.off()
  
}

barplot(dpois(seq(0,10,1),lambda=0.000625))

#analyze the fsts--as calculated by C++
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
pcr.files<-list.files(pattern="summstatspcr")

png("fsts_summstats_pcr.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(6,1),oma=c(2,2,2,2),mar=c(2,2,2,2))
for(i in 1:length(pcr.files)){
  summ<-data.frame(read.table(pcr.files[i],sep='\t',
  	header=T,stringsAsFactors=F))
  summ$Fst<-as.numeric(summ$Fst)
  ifsts<-summ[!is.na(summ$Fst),]
  ifsts<-ifsts[as.numeric(ifsts$Fst)>=0,]
  
  scaffs<-levels(as.factor(ifsts[,1]))
  scaffs[1:22]<-lgs
  
  i1<-fst.plot(ifsts,ci.dat=c(-10,10),sig.col=c("black","black"),fst.name="Fst",
  	chrom.name="Chrom",bp.name="Pos",axis.size=1,
  	groups=as.factor(scaffs[scaffs %in% levels(factor(ifsts[,1]))]))
  last<-0
  for(ii in 1:length(lgs)){
    text(x=mean(i1[i1$Chrom ==lgs[ii],"Pos"]),y=-0.01,
         labels=lgs[ii], adj=1, xpd=TRUE,srt=90,cex=1)
    last<-max(i1[i1$Chrom ==lgs[ii],"Pos"])
  }
  mtext(pcr.files[i],3)
}
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

#analyze the fsts--using fst.two.vcf
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
vcf.files<-list.files(pattern="pcr.\\d+.vcf")

png("fsts_vcf_pcr.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(6,1),oma=c(2,2,2,2),mar=c(2,2,2,2))
for(i in 1:length(vcf.files)){
  vcf<-parse.vcf(vcf.files[i])
  vcf$SNP<-paste(vcf$`#CHROM`, vcf$POS,sep=".")
  sd.ind<-colnames(vcf)[grep("sd",colnames(vcf))]
  dd.ind<-colnames(vcf)[grep("dd",colnames(vcf))]
  locus.info<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
  sd.vcf<-cbind(vcf[,locus.info],vcf[,sd.ind])
  dd.vcf<-cbind(vcf[,locus.info],vcf[,dd.ind])
  fst<-do.call("rbind",apply(sd.vcf,1,fst.two.vcf,vcf2=dd.vcf,match.index="SNP"))
  keep.fst<-fst[fst$Num1 > 0 & fst$Num2 > 0,]
  scaffs<-levels(as.factor(keep.fst[,1]))
  scaffs[1:22]<-lgs
  
  v1<-fst.plot(keep.fst,ci.dat=c(-10,10),sig.col=c("black","black"),fst.name="Fst",
               chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% levels(factor(keep.fst[,1]))]))
  last<-0
  for(ii in 1:length(lgs)){
    text(x=mean(v1[v1$Chrom ==lgs[ii],"Pos"]),y=-0.01,
         labels=lgs[ii], adj=1, xpd=TRUE,srt=90,cex=1)
    last<-max(v1[v1$Chrom ==lgs[ii],"Pos"])
  }
  mtext(vcf.files[i],3,cex=0.75)
}
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

#Effect of Ne (migration rates)
ne.files<-c(list.files(pattern="summstats.ne\\d+"),list.files(pattern="summstats.u"))

png("fsts_mutation.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(6,1),oma=c(2,2,2,2),mar=c(2,2,2,2))
for(i in 1:length(ne.files)){
    summ<-data.frame(read.table(ne.files[i],sep='\t',
                                header=T,stringsAsFactors=F))
    summ$Fst<-as.numeric(summ$Fst)
    ifsts<-summ[!is.na(summ$Fst),]
    ifsts<-ifsts[as.numeric(ifsts$Fst)>=0,]
    
    scaffs<-levels(as.factor(ifsts[,1]))
    scaffs[1:22]<-lgs
    
    i1<-fst.plot(ifsts,ci.dat=c(-10,10),sig.col=c("black","black"),fst.name="Fst",
                 chrom.name="Chrom",bp.name="Pos",axis.size=1,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(ifsts[,1]))]))
    last<-0
    for(ii in 1:length(lgs)){
      text(x=mean(i1[i1$Chrom ==lgs[ii],"Pos"]),y=-0.01,
           labels=lgs[ii], adj=1, xpd=TRUE,srt=90,cex=1)
      last<-max(i1[i1$Chrom ==lgs[ii],"Pos"])
    }
    mtext(pcr.files[i],3)
  }
  mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

afs.vcf<-parse.vcf("ssc_insilico.afs.vcf")
afs.afs<-do.call("rbind",apply(afs.vcf,1,calc.afs.vcf))

afs.list<-list.files(pattern="summstats.afs")
