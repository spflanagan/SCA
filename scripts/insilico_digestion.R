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

summstats.plot<-function(file.list,make.png=T,file.name=""){
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
pcr.files<-list.files(pattern="summstatspcr")
summstats.plot(pcr.files,file.name="fsts_summstats_pcr.png")

#Effect of Ne (migration rates)
ne.files<-c(list.files(pattern="summstats.ne\\d+"),list.files(pattern="summstats.u"))
summstats.plot(ne.files,file.name="fsts_mutation.png")

#skewed allele frequency
afs.list<-list.files(pattern="summstats.afs.pcr\\d+.txt")
summstats.plot(afs.list,file.name="fsts_afs.png")

noshear<-read.table("ssc_insilico_summstats.afs.pcr02.noshear.txt",sep='\t',
                    header=T,stringsAsFactors=F)

null<-read.table("ssc_insilico_summstats.null.skewedAFS.txt",sep='\t',
                 header=T,stringsAsFactors=F)
##Rearrange data for plots
pcr.afs.dat<-data.frame(PCR=numeric(),AFS=numeric(),Fst=numeric(),stringsAsFactors = FALSE)
for(i in 1:(length(pcr.files)-1)){
  dat<-data.frame(read.delim(pcr.files[i],stringsAsFactors=FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1])
  fsts<-fsts[!is.na(fsts)]
  r<-gsub("ssc_insilico_summstatspcr.(\\d+).txt","\\1",pcr.files[i])
  pcr.afs.dat<-data.frame(rbind(pcr.afs.dat,
    cbind(PCR=rep(as.numeric(r),length(fsts)),AFS=rep(as.numeric(0.5),length(fsts)),Fst=fsts)),stringsAsFactors=F)
}
for(i in 1:length(afs.list)){
  dat<-data.frame(read.delim(afs.list[i],stringsAsFactors=FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1])
  fsts<-fsts[!is.na(fsts)]
  r<-gsub("ssc_insilico_summstats.afs.pcr.(\\d+).txt","\\1",afs.list[i])
  pcr.afs.dat<-data.frame(rbind(pcr.afs.dat,
    cbind(PCR=rep(as.numeric(r),length(fsts)),AFS=rep(as.numeric(0.8),length(fsts)),Fst=fsts)),stringsAsFactors=F)
}

mut.dat<-data.frame(Ne=numeric(),MinMut=character(),Fst=numeric(),stringsAsFactors=FALSE)
for(i in 1:length(ne.files))
{
  dat<-data.frame(read.delim(ne.files[i],stringsAsFactors = FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1])
  fsts<-fsts[!is.na(fsts)]
  if(length(grep("ne",ne.files[i]))>0){
    n<-gsub("ssc_insilico_summstats.ne(\\d+).txt","\\1",ne.files[i])
    m<-"10E8"
  }
  if(length(grep("Ne",ne.files[i]))>0){ 
    n<-gsub("ssc_insilico_summstats.\\w+Ne(\\d+).txt","\\1",ne.files[i])
    m<-"10E7"
  }
  
  mut.dat<-data.frame(rbind(mut.dat,
    cbind(Ne=rep(as.numeric(n),length(fsts)),MinMut=rep(as.character(m),length(fsts)),Fst=as.numeric(fsts))),stringsAsFactors=FALSE)
}
mut.dat$Fst<-as.numeric(as.character(mut.dat$Fst))

###PLOTS####
boxplot(pcr.afs.dat$Fst~pcr.afs.dat$AFS*pcr.afs.dat$PCR,
  col=c("white","grey"),pch=19,axes=F,ylim=c(-0.02,0.06))
axis(1,at=c(0,seq(1.5,11.5,2)),labels=c("",seq(1,5),""),pos=-0.02)
axis(2,pos=0.1,las=1,at=seq(-0.02,0.06,0.01))
legend("topleft",c("Mean p = 0.5","Mean p = 0.8"),pt.bg=c("white","grey"),
       pch=22,bty="n",ncol=2)

boxplot(mut.dat$Fst~mut.dat$MinMut*as.numeric(as.character(mut.dat$Ne)),
  col=c("white","grey"),pch=21,axes=F,ylim=c(-0.06,0.06),
  outcol="black",outbg=c("white","grey"))
axis(1,at=c(0,seq(1.5,8.5,2)),labels=c("",5000,10000,20000,""),pos=-0.06)
axis(2,pos=0.3,las=1,at=seq(-0.06,0.06,0.01))
legend("topleft",
  c(expression(mu~"from ["~10^{-9}~","~10^{-8}~"]"),expression(mu~"from ["~10^{-8}~","~10^{-7}~"]")),
  pt.bg=c("white","grey"),
  pch=21,bty="n")

hist(log(null$Fst))

####FOR COMPARISON####
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

