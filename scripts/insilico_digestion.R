#Author: Sarah P. Flanagan
#Last Updated: 28 September 2016
#Date Started: 23 September 2016
#Purpose: analyze in-silico digestion results

rm(list=ls())
setwd("~/Projects/SCA/results/insilico")
#setwd("//VBOXSVR/Home/Projects/SCA/results/insilico")
source("../../scripts/plotting_functions.R")

#sd<-read.delim("SSC_sdigested.txt",header=T)
#sd$RS<-paste(sd$Chrom,sd$PstIPos,sep=".")
#dd<-read.delim("SSC_ddigested.txt",header=T)
#dd$RS<-paste(dd$Chrom,dd$PstIPos,sep=".")
#hist(sd$FragmentLength)

summstats.plot<-function(file.list,make.png=T,file.name=""){
  lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
         "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
         "LG22")
  if(make.png==T) png(file.name,height=10,width=7.5,units="in",res=300)
  par(mfrow=c(length(file.list),1),oma=c(2,2,2,2),mar=c(2,2,2,2))
  for(i in 1:length(file.list)){
    summ<-data.frame(read.table(file.list[i],sep='\t',
                                header=T,stringsAsFactors=F))
    summ$Fst[summ$Fst=="-nan(ind)"]<-NA
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

#barplot(dpois(seq(0,10,1),lambda=0.000625))

#analyze the fsts--as calculated by C++
#list files
null<-read.delim("ssc_insilico_summstats.null.txt",stringsAsFactors = F)
null.files<-list.files(pattern="summstats.null")
pcr.files<-list.files(pattern="summstats\\.pcr")
mut.files<-list.files(pattern="summstats.u")
shear.files<-list.files(pattern="summstats.shearingbias")
comb.files<-list.files(pattern="summstats.comb")

#Summstats plots
summstats.plot(pcr.files,file.name="fsts_summstats_pcr.png")
summstats.plot(mut.files,file.name="fsts_mutation.png")
summstats.plot(shear.files,file.name="fsts_shear.png")


##Rearrange data for plots
#for panel 1
pcr.afs.dat<-data.frame(PCR=numeric(),AFS=numeric(),Fst=numeric(),stringsAsFactors = FALSE)
for(i in 1:length(pcr.files)){
  dat<-data.frame(read.delim(pcr.files[i],stringsAsFactors=FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1 & dat$Fst!="-nan(ind)"])
  fsts<-fsts[!is.na(fsts)]
  r<-gsub("ssc_insilico_summstats.pcr(\\d+).\\w+.txt","\\1",pcr.files[i])
  s<-gsub("ssc_insilico_summstats.pcr\\d+.(\\w+).txt","\\1",pcr.files[i])
  if(s == "notskewed") afs<-0.5
  if(s == "skewed") afs <- 0.8
  pcr.afs.dat<-data.frame(rbind(pcr.afs.dat,
    cbind(PCR=rep(as.numeric(r),length(fsts)),AFS=rep(as.numeric(afs),length(fsts)),Fst=fsts)),stringsAsFactors=F)
}
#add null info
dat<-data.frame(read.delim("ssc_insilico_summstats.null.skewedAFS.txt",stringsAsFactors=FALSE))
fsts<-as.numeric(dat$Fst[dat$Fst!=-1 & dat$Fst!="-nan(ind)"])
fsts<-fsts[!is.na(fsts)]
pcr.afs.dat<-data.frame(rbind(pcr.afs.dat,
    cbind(PCR=rep(as.numeric(0),length(fsts)),AFS=rep(as.numeric(0.8),length(fsts)),Fst=fsts)),stringsAsFactors=F)
fsts<-as.numeric(null$Fst[null$Fst!=-1])
fsts<-fsts[!is.na(fsts)]
pcr.afs.dat<-data.frame(rbind(pcr.afs.dat,
  cbind(PCR=rep(as.numeric(0),length(fsts)),AFS=rep(as.numeric(0.5),length(fsts)),Fst=fsts)),stringsAsFactors=F)

#for panel 2
mut.dat<-data.frame(Ne=numeric(),MinMut=character(),Fst=numeric(),stringsAsFactors=FALSE)
for(i in 1:length(mut.files))
{
  dat<-data.frame(read.delim(mut.files[i],stringsAsFactors = FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1 & dat$Fst!="-nan(ind)"])
  fsts<-fsts[!is.na(fsts)]
  n<-gsub("ssc_insilico_summstats.u10E\\dNe(\\d+).txt","\\1",mut.files[i])
  m<-gsub("ssc_insilico_summstats.(u10E\\d)Ne\\d+.txt","\\1",mut.files[i])
  mut.dat<-data.frame(rbind(mut.dat,
    cbind(Ne=rep(as.numeric(n),length(fsts)),MinMut=rep(as.character(m),length(fsts)),Fst=as.numeric(fsts))),stringsAsFactors=FALSE)
}
mut.dat$Fst<-as.numeric(as.character(mut.dat$Fst))
#for panel 3
p3.files<-c(null.files,shear.files)
null.dat<-data.frame(Shearing=character(),Asymmetric=numeric(),Skew=character(),Fst=numeric(),stringsAsFactors=FALSE)
for(i in 1:length(p3.files))
{
  dat<-data.frame(read.delim(p3.files[i],stringsAsFactors = FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1 & dat$Fst!="-nan(ind)"])
  fsts<-fsts[!is.na(fsts)]
  if(length(grep("shearingbias",p3.files[i]))>0){
    a<-"ShearingBias"
  }else{
    a<-"NoShearingBias"
  }
  if(length(grep("asymmetric",p3.files[i]))>0){
    b<-"Asymmetric"
  }else{
    b<-"Symmetric"
  }
  if(length(grep("notskewed",p3.files[i]))>0){ 
    c<-0.5
  }else{
    if(length(grep("skewed",p3.files[i]))>0){
      c<-0.8}
    else{
      c<-0.5}
  }
  print(paste(i,a,b,c))
  null.dat<-data.frame(rbind(null.dat,
    cbind(Shearing=rep(as.character(a),length(fsts)),
          Asymmetric=rep(as.character(b),length(fsts)),Skew=rep(as.numeric(c),length(fsts)),Fst=as.numeric(fsts))),stringsAsFactors=FALSE)
}
null.dat$Fst<-as.numeric(as.character(null.dat$Fst))

comb.dat<-data.frame(Description=character(),Fst=numeric(),stringsAsFactors=FALSE)
for(i in 1:length(comb.files))
{
  dat<-data.frame(read.delim(comb.files[i],stringsAsFactors = FALSE))
  fsts<-as.numeric(dat$Fst[dat$Fst!=-1 & dat$Fst!="-nan(ind)"])
  fsts<-fsts[!is.na(fsts)]
  name<-gsub("ssc_insilico_summstats.comb.(\\w+.*).txt","\\1",comb.files[i])
  comb.dat<-rbind(comb.dat,cbind(Description=rep(name,length(fsts)),
    Fst=fsts))
}
comb.dat$Fst<-as.numeric(as.character(comb.dat$Fst))
###PLOTS####
jpeg("insilico_results_revised.jpeg",height=8,width=10.5,units="in",res=300)
par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
boxplot(pcr.afs.dat$Fst~pcr.afs.dat$AFS*pcr.afs.dat$PCR,
  col=c("white","grey"),pch=21,axes=F,ylim=c(0,0.02),
  outcol="black",outbg=c("white","grey"))
axis(1,at=c(0,seq(0.5,12.5,2)),labels=rep("",8),pos=-0.0008,tick=T)
axis(1,at=seq(1.5,11.5,2),labels=seq(0,5),pos=-0.00001,tck=0,tick=F)
axis(2,pos=0.5,las=1,at=seq(-0.02,0.06,0.01))
legend("top",c("Mean p = 0.5","Mean p = 0.8"),pt.bg=c("white","grey"),
       pch=21,bty="n",ncol=2)
mtext("% of reads due to PCR duplicates per cycle",1,outer=F,line=1.5,cex=0.75)

boxplot(mut.dat$Fst~mut.dat$MinMut*as.numeric(as.character(mut.dat$Ne)),
  col=c("white","grey"),pch=21,axes=F,ylim=c(-0.2,0.2),
  outcol="black",outbg=c("white","grey"))
axis(1,at=seq(-1.5,8.5,2),labels=rep("",6),pos=-0.2)
axis(1,at=seq(1.5,5.5,2),labels=c(5000,10000,20000),pos=-0.18,tick=F,tck=0)
axis(2,pos=0.5,las=1,at=seq(-0.5,0.3,0.1))
legend("top",ncol=2,
  c(expression(mu~"from ["~10^{-8}~","~10^{-7}~"]"),expression(mu~"from ["~10^{-9}~","~10^{-8}~"]")),
  pt.bg=c("white","grey"),
  pch=21,bty="n")
mtext(expression(italic(N)[italic(e)]),1,outer=F,line=1,cex=0.75)

boxplot(null.dat$Fst~factor(null.dat$Asymmetric)*factor(null.dat$Skew)*factor(null.dat$Shearing),
        col=c("white","grey"),pch=c(rep(21,4),rep(23,4)),axes=F,ylim=c(0,0.03),
        outcol="black",outbg=c("white","grey"))
axis(1,at=seq(-1.5,8.5,2),labels=rep("",6),pos=-0.001)
axis(1,at=seq(1.5,7.5,2),labels=rep(c("Mean p = 0.5","Mean p = 0.8"),2),pos=-0.00025,tick=F,tck=0)
axis(1,at=c(2.5,6.5),labels=c("No Shearing Bias","Shearing Bias"),pos=-0.002,tick=F,tck=0)
axis(2,pos=0.5,las=1,at=seq(-0.01,0.03,0.01))
legend("top",c("Asymmetric Sampling","Symmetric Sampling","No Shearing Bias","Shearing Bias"),
       pt.bg=c("white","grey","white","white"),
       pch=c(22,22,21,23),bty="n",ncol=2)

boxplot(comb.dat$Fst~comb.dat$Description,axes=F,col=alpha(rainbow(5),0.5),
        outbg=alpha(rainbow(5),0.5),pch=21,ylim=c(-0.04,0.04))
axis(1,at=c(-0.5,0.5,5.5,10.5,15.5,17),labels=rep("",6),pos=-0.04)
axis(1,at=c(3,8,13),pos=-0.04,tick=F,tck=0,
     labels=c("1% of reads from\nPCR duplicates","2% of reads from\nPCR duplicates","3% of reads from\nPCR duplicates"))
axis(2,pos=0.5,las=1,at=seq(-0.04,0.04,0.02))
legend("top",c(expression(mu~"from ["~10^{-8}~","~10^{-7}~"]"),"Mean p = 0.5",
                       expression(mu~"from ["~10^{-9}~","~10^{-8}~"]"),"No Shearing Bias","Symmetric Sampling"),
       pt.bg=c(alpha(rainbow(5),0.5)),pch=21,bty="n",ncol=2)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

##STATS##
summary(aov(pcr.afs.dat$Fst~pcr.afs.dat$AFS*pcr.afs.dat$PCR))
aov(mut.dat$Fst~mut.dat$MinMut*as.numeric(as.character(mut.dat$Ne)))
aov(null.dat$Fst~factor(null.dat$Asymmetric)*factor(null.dat$Skew)*factor(null.dat$Shearing))

comb.dat$PCR<-gsub("pcr(\\d).\\w+.*","\\1",comb.dat$Description)
comb.dat$MutationRate<-gsub("\\w+.*(u10E\\dNe10000).\\w+.*","\\1",comb.dat$Description)
comb.dat$Symmetric<-gsub("\\w+.*10000.(\\w?symmetric).\\w+.*","\\1",comb.dat$Description)
comb.dat$Shearing<-gsub("\\w+.*symmetric\\.(\\w+\\.){0,}(\\w+)","\\2",comb.dat$Description)
comb.dat$Shearing[comb.dat$Shearing=="skewed"]<-"nobias"
comb.dat$Skew<-gsub("\\w+.*symmetric\\.(\\w+)","\\1",comb.dat$Description)
comb.dat$Skew[comb.dat$Skew=="skewed.shearingbias"]<-"skewed"
comb.dat$Skew[comb.dat$Skew=="shearingbias"]<-"notskewed"
comb.dat$PCR<-as.factor(comb.dat$PCR)
comb.dat$MutationRate<-as.factor(comb.dat$MutationRate)
comb.dat$Symmetric<-as.factor(comb.dat$Symmetric)
comb.dat$Shearing<-as.factor(comb.dat$Shearing)
comb.dat$Skew<-as.factor(comb.dat$Skew)
summary(lm(Fst~PCR*MutationRate*Symmetric+Shearing+Skew,data = comb.dat))

####FOR COMPARISON####
#analyze the fsts--using fst.two.vcf
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
vcf.files<-list.files(pattern="u10E\\dNe\\d+.vcf")

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

