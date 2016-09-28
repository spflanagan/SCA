#Author: Sarah P. Flanagan
#Last updated: 23 September 2016
#Date Started: 24 February 2016
#Purpose: compare RAD libraries prepared by two different methods

rm(list=ls())
library(ggplot2)
setwd("~/Projects/SCA/results")
source("../scripts/plotting_functions.R")

###############FUNCTIONS####################
parse.vcf<-function(filename){
  vcf<-read.delim(filename,comment.char="#",sep='\t',header=F,stringsAsFactors = F)
  header.start<-grep("#CHROM",scan(filename,what="character"))
  header<-scan(filename,what="character")[header.start:(header.start+ncol(vcf)-1)]
  colnames(vcf)<-header
  return(vcf)
}

vcf.cov.loc<-function(vcf.row,subset){
  cov<-unlist(lapply(vcf.row[subset],function(x){ 
    c<-strsplit(as.character(x),split=":")[[1]][3]
    return(c)
  }))
  miss<-length(cov[cov==".,."])
  pres<-length(cov[cov!=".,."])
  ref<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
    function(x){
      strsplit(as.character(x),",")[[1]][1] 
    }))))/pres
  alt<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
    function(x){
      strsplit(as.character(x),",")[[1]][2] 
    }))))/pres
  tot<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
    function(x){
      as.numeric(strsplit(as.character(x),",")[[1]][1]) + 
      as.numeric(strsplit(as.character(x),",")[[1]][2])
    }))))
  var.cov<-var(as.numeric(unlist(lapply(cov[cov!=".,."],
                                    function(x){
                                      as.numeric(strsplit(as.character(x),",")[[1]][1]) + 
                                        as.numeric(strsplit(as.character(x),",")[[1]][2])
                                    }))))
  het<-unlist(lapply(vcf.row[subset],function(x){ 
    strsplit(as.character(x),split=":")[[1]][1]
  }))
  het<-length(het[het=="0/1" | het=="1/0"])
  return(data.frame(Chrom=vcf.row[1],Pos=vcf.row["POS"],Locus=vcf.row["ID"],
    NumMissing=miss, NumPresent=pres,PropMissing=miss/(miss+pres),
    AvgCovRef=ref,AvgCovAlt=alt, AvgCovRatio=ref/alt,AvgCovTotal=tot/pres, CovVariance=var.cov,
    NumHet=het,PropHet=het/pres,TotalNumReads = tot,stringsAsFactors = F))
}

vcf.cov.ind<-function(vcf.col){
  cov<-unlist(lapply(vcf.col,function(x){ 
    c<-strsplit(as.character(x),split=":")[[1]][3]
    return(c)
  }))
  miss<-length(cov[cov==".,."])
  pres<-length(cov[cov!=".,."])
  ref<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                    function(x){
                                      strsplit(as.character(x),",")[[1]][1] 
                                    }))))/pres
  alt<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                    function(x){
                                      strsplit(as.character(x),",")[[1]][2] 
                                    }))))/pres
  tot<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                    function(x){
                                      as.numeric(strsplit(as.character(x),",")[[1]][1]) + 
                                        as.numeric(strsplit(as.character(x),",")[[1]][2])
                                    }))))/pres
  het<-unlist(lapply(vcf.col,function(x){ 
    strsplit(as.character(x),split=":")[[1]][1]
  }))
  het<-length(het[het=="0/1" | het=="1/0"])
  return(list(NumMissing=miss,NumPresent=pres,AvgCovRef=ref,AvgCovAlt=alt,AvgCovTot=tot,PropHet=het/pres, NumReads=tot*pres))
}
  
fst.two.vcf<-function(vcf1.row,vcf2,match.index){
  #match.index is the column used to match the two
  #use in conjunction with apply
    #e.g. apply(vcf,1,fst.two.vcf,vcf2=vcf.2,match.index="SNP")
  hs1<-hs2<-hs<-ht<-0
  freqall<-gt1<-gt2<-NULL
  vcf2.row<-vcf2[vcf2[,match.index]%in%vcf1.row[match.index],]
  if(nrow(vcf2.row)>1)#first make sure we have one reading per locus
  {
    print("Multiple instances in vcf2.")
    fst<-NA
  }
  else{
    if(nrow(vcf2.row)==0)
    {
      print("No instances in vcf2.")
      fst<-NA
    }else #we're good to go
    {
      gt1<-unlist(lapply(vcf1.row,function(x){ 
        c<-strsplit(as.character(x),split=":")[[1]][1]
        return(c)
      }))
      gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
      gt1[gt1=="1/0"]<-"0/1"
      gt1<-gsub(pattern = "0",replacement = vcf1.row["REF"],gt1)
      gt1<-gsub(pattern = "1",replacement = vcf1.row["ALT"],gt1)
      al1<-unlist(strsplit(as.character(gt1),split = "/"))
      gt2<-unlist(lapply(vcf2.row,function(x){ 
        c<-strsplit(as.character(x),split=":")[[1]][1]
        return(c)
      }))
      gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
      gt2[gt2=="1/0"]<-"0/1"
      gt2<-gsub(pattern = "0",replacement = vcf2.row["REF"],gt2)
      gt2<-gsub(pattern = "1",replacement = vcf2.row["ALT"],gt2)
      al2<-unlist(strsplit(as.character(gt2),split="/"))
       #calculate frequencies
      freq1<-summary(factor(al1))/sum(summary(factor(al1)))	
      freq2<-summary(factor(al2))/sum(summary(factor(al2)))	
      freqall<-summary(as.factor(c(al1,al2)))/
        sum(summary(as.factor(c(al1,al2))))
      hets<-c(names(freq1)[2],names(freq2)[2])
      if(length(freq1)>1){ hs1<-1-sum(freq1^2) 
      } else {
        hs1<-0
      }
      if(length(freq2)>1){ hs2<-1-sum(freq2^2)
      } else {
        hs2<-0
      }
      if(length(freqall)>1){
        hs<-mean(c(hs1,hs2))
        ht<-1-sum(freqall^2)
        fst<-(ht-hs)/ht
      }
      if(length(freqall)<=1){ fst<-1 }
    }#end else good to go
  }#end else vcf2

  return(data.frame(Chrom=vcf1.row["#CHROM"],Pos=vcf1.row["POS"],
    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
    Num1=length(gt1),Num2=(length(gt2))))
}#end function

######FILES#####
drad<-parse.vcf("drad.vcf")
orad<-parse.vcf("orad.vcf")
both<-parse.vcf("both.vcf")

o.ind<-grep("orad",colnames(both),value=T)
d.ind<-grep("sample",colnames(both),value=T)
#########ANALYSIS##########
###coverage per locus
#assembled together

bo.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=o.ind))
bd.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=d.ind))
b.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=c(o.ind,d.ind)))
loc.cov<-merge(o.cov,d.cov,by="Locus")
#because some of them have no alt
loc.cov$oRAD.AvgCovRatio[loc.cov$oRAD.AvgCovRatio=="Inf"]<-loc.cov[loc.cov$oRAD.AvgCovRatio=="Inf","AvgCovRef.x"]
loc.cov$dRAD.AvgCovRatio[loc.cov$dRAD.AvgCovRatio=="Inf"]<-loc.cov[loc.cov$dRAD.AvgCovRatio=="Inf","AvgCovRef.y"]

t.test(loc.cov$oRAD.AvgCovRatio,loc.cov$dRAD.AvgCovRatio,paired=T,alternative="less")
#private alleles
dim(loc.cov[(loc.cov$AvgCovRef.x == 0 & loc.cov$AvgCovAlt.y == 0) | (loc.cov$AvgCovAlt.x == 0 & loc.cov$AvgCovRef.y == 0),])#2042
dim(loc.cov[loc.cov$NumHet.x==0,])#72307
dim(loc.cov[loc.cov$NumHet.y==0,])#14436

#assembled separately
o.cov<-do.call("rbind",apply(orad,1,vcf.cov.loc,subset=o.ind))
d.cov<-do.call("rbind",apply(drad,1,vcf.cov.loc,subset=d.ind))
#compare assembly methods
wilcox.test(loc.cov$oRAD.AvgCovRatio,o.cov$AvgCovRatio,"greater")
wilcox.test(loc.cov$dRAD.AvgCovRatio,d.cov$AvgCovRatio,"less")
wilcox.test(loc.cov$AvgCovTotal.x,o.cov$AvgCovTotal,"less")
wilcox.test(loc.cov$AvgCovTotal.y,d.cov$AvgCovTotal,"greater")
#compare sd vs dd
lc.comp<-data.frame(LibraryPrep=c(rep("sdRAD",nrow(o.cov)),rep("ddRAD",nrow(d.cov)),
    rep("sdRAD",nrow(loc.cov)),rep("ddRAD",nrow(loc.cov))), 
  Assembly=c(rep("Alone",nrow(o.cov)),rep("Alone",nrow(d.cov)),
             rep("Together",nrow(loc.cov)),rep("Together",nrow(loc.cov))),
  AvgCovRatio=c(o.cov$AvgCovRatio,d.cov$AvgCovRatio,loc.cov$oRAD.AvgCovRatio,loc.cov$dRAD.AvgCovRatio),
  AvgCovTotal=c(o.cov$AvgCovTotal,d.cov$AvgCovTotal,loc.cov$AvgCovTotal.x,loc.cov$AvgCovTotal.y))
#NumReadsPerInd=c(oc$NumReadsAlone,dc$NumReadsAlone,oc$NumReadsTogether,dc$NumReadsTogether)
#using separate assemblies
wilcox.test(o.cov$AvgCovRatio,d.cov$AvgCovRatio,"less")
wilcox.test(o.cov$AvgCovTotal,d.cov$AvgCovTotal)
#using same assembly
wilcox.test(loc.cov$oRAD.AvgCovRatio,loc.cov$dRAD.AvgCovRatio,paired=T,"less")
wilcox.test(loc.cov$AvgCovTotal.x,loc.cov$AvgCovTotal.y,paired=T,"less")
summary(aov(lc.comp$AvgCovRatio~lc.comp$LibraryPrep*lc.comp$Assembly))

jpeg("CoverageComparisonBoxplot.jpeg",height=7,width=7,units="in",res=300)
par(mfrow=c(2,1),oma=c(1,2,1,1),mar=c(2,2,2,2))
boxplot(log(lc.comp$AvgCovRatio+1)~lc.comp$LibraryPrep*lc.comp$Assembly,names=F,
        xaxt='n',col=c("cadetblue","coral1"))
axis(1,at=c(1.5,3.5),labels=c("Alone","Together"))
mtext("log(Avg Ratio Num Reads\nRef/Alt)",2,outer=F,line=2)
legend("topleft",c("sdRAD","ddRAD"),col=c("cadetblue","coral1"),pch=15,bty='n')
boxplot(log(lc.comp$AvgCovTotal+1)~lc.comp$LibraryPrep*lc.comp$Assembly,names=F,
        xaxt='n',col=c("cadetblue","coral1"))
mtext("log(Avg Num Reads\nPer Individual Per Locus)",2,outer=F,line=2)
axis(1,at=c(1.5,3.5),labels=c("Alone","Together"))
dev.off()

#Coverage by individual--assembled in one
o.icov<-as.data.frame(do.call("rbind",apply(both[,o.ind],2,vcf.cov.ind)))
o.icov$Method<-rep("oRAD",nrow(o.icov))
d.icov<-as.data.frame(do.call("rbind",apply(both[,d.ind],2,vcf.cov.ind)))
d.icov$Method<-rep("dRAD",nrow(d.icov))
ind.cov<-as.data.frame(rbind(o.icov,d.icov))
ind.cov$NumReads<-as.numeric(ind.cov$AvgCovTot)*as.numeric(ind.cov$NumPresent)
ic<-as.data.frame(apply(ind.cov,2,as.numeric))
rownames(ic)<-rownames(ind.cov)
ic$Method<-as.factor(ind.cov$Method)
#Coverage by individual--assembled separately
orad.icov<-as.data.frame(do.call("rbind",apply(orad[,10:ncol(orad)],2,vcf.cov.ind)))
oic<-apply(orad.icov,2,as.numeric)
rownames(oic)<-rownames(orad.icov)
drad.icov<-as.data.frame(do.call("rbind",apply(drad[,10:ncol(drad)],2,vcf.cov.ind)))
dic<-apply(drad.icov,2,as.numeric)
rownames(dic)<-rownames(drad.icov)

#write these to file so I can easily pick back up later
write.table(oic,"orad_coverage.csv",quote=T,row.names=T,col.names=T,sep='\t')
write.table(dic,"drad_coverage.csv",quote=T,row.names=T,col.names=T,sep='\t')
write.table(ic,"both_coverage.csv",quote=T,row.names=T,col.name=T,sep='\t')

oic<-read.delim("orad_coverage.csv")
dic<-read.delim("drad_coverage.csv")
ic<-read.delim("both_coverage.csv")

ic.o<-ic[rownames(ic) %in% rownames(oic),]
ic.d<-ic[rownames(ic) %in% rownames(dic),]

mo<-merge(ic.o,oic,by=0)
md<-merge(ic.d,dic,by=0)
oc<-data.frame(NumReadsTogether=mo$NumReads.x,NumReadsAlone=mo$NumReads.y,row.names=mo$Row.names)
dc<-data.frame(NumReadsTogether=md$NumReads.x,NumReadsAlone=md$NumReads.y,row.names=md$Row.names)
t.test(oc$NumReadsTogether,oc$NumReadsAlone,paired=T,alternative="less")
t.test(dc$NumReadsTogether,dc$NumReadsAlone,paired=T,alternative="greater")

#plot
jpeg("CoverageAssemblyMethodComp.jpeg",height=10.5,width=7,units="in",res=300)
par(mfrow=c(3,2),oma=c(1,1,1,1),mar=c(2,2,2,2))
hist(log(loc.cov$AvgCovTotal.x),col=rgb(0,0,1,0.5),main="sdRAD-seq",axes=F,xlab="",ylab="",
     xlim=c(0,10.5),ylim=c(0,200000),breaks=seq(0,10.5,0.5))
hist(log(o.cov$AvgCovTotal), col=rgb(1,0,0,0.5), add=T,breaks=seq(0,10.5,0.5))
axis(1,pos=0)
axis(2,pos=0)
mtext("Number of Loci",2,outer=F,line=1,cex=0.75)
hist(log(loc.cov$AvgCovTotal.y), col=rgb(0,0,1,0.5),main="ddRAD-seq",axes=F,ylab="",xlab="",
     xlim=c(0,8),ylim=c(0,150000),breaks=seq(0,8,0.5))
hist(log(d.cov$AvgCovTotal), col=rgb(1,0,0,0.5), add=T,breaks=seq(0,8,0.5))
axis(1,pos=0)
axis(2,pos=0)
mtext("ln(Average Number of Reads Per Individual)",1,outer=T,line=-52,cex=0.75)


hist(log(loc.cov$oRAD.AvgCovRatio),col=rgb(0,0,1,0.5),main="",axes=F,xlab="",ylab="",breaks=seq(-10,10,1),ylim=c(0,150000))
hist(log(o.cov$AvgCovRatio), col=rgb(1,0,0,0.5), add=T,breaks=seq(-10,10,1))
axis(1,pos=0)
axis(2,pos=-10)
mtext("Number of Loci",2,outer=F,line=1,cex=0.75)
hist(log(loc.cov$dRAD.AvgCovRatio), col=rgb(0,0,1,0.5),main="",axes=F,ylab="",xlab="",breaks=seq(-5,14,1),ylim=c(0,150000))
hist(log(d.cov$AvgCovRatio), col=rgb(1,0,0,0.5), add=T,breaks=seq(-5,14,1))
axis(1,pos=0)
axis(2,pos=-5)
mtext("ln(Average Number of Reads in Ref/Avg Number of Reads in Alt)",1,outer=T,line=-27,cex=0.75)


hist(oc$NumReadsTogether, col=rgb(0,0,1,0.5), ylim=c(0,30), breaks=seq(200000,8000000,400000),
     main="",axes=F,xlim=c(0,8000000))
hist(oc$NumReadsAlone,col=rgb(1,0,0,0.5), add=T,breaks=seq(200000,8000000,400000))
axis(1,pos=0)
axis(2,pos=0)
mtext("Number of Individuals",2,outer=F,line=1,cex=0.75)
hist(dc$NumReadsTogether, col=rgb(0,0,1,0.5), ylim=c(0,150), breaks=seq(2000,4000000,200000),
     main="",axes=F)
hist(dc$NumReadsAlone,col=rgb(1,0,0,0.5), add=T,breaks=seq(2000,4000000,200000))
axis(1,pos=0)
axis(2,pos=0)
legend("topright",c("Assembled Together", "Assembled Separately"),pch=15,col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),bty='n')
mtext("Total Number of Reads",1,outer=T,cex=0.75)
dev.off()

wilcox.test(oc$NumReadsTogether,oc$NumReadsAlone,paired=T,alternative="less")
wilcox.test(dc$NumReadsTogether,dc$NumReadsAlone,paired=T,alternative="greater")


###Variance in coverage
#assembled separately
hist(log(o.cov$CovVariance),col=alpha("cadetblue",0.5))
hist(log(d.cov$CovVariance),col=alpha("coral1",0.5),add=T)
boxplot(log(o.cov$CovVariance+1),log(d.cov$CovVariance+1),col=c("cadetblue","coral1"),
        names=c("sdRAD","ddRAD"),ylab="Variance in Coverage",ylim=c(0,22))
text(x=c(1,2),y=c(21,21),c("250425","69109"))
wilcox.test(log(o.cov$CovVariance+1),log(d.cov$CovVariance+1),"greater")

d.sub<-d.cov[sample(nrow(d.cov),30000,replace=F),]
o.sub<-o.cov[sample(nrow(o.cov),30000,replace=F),]
wilcox.test(log(o.sub$CovVariance+1),log(d.sub$CovVariance+1),"greater")
boxplot(log(o.sub$CovVariance+1),log(d.sub$CovVariance+1),col=c("cadetblue","coral1"),
        names=c("sdRAD","ddRAD"),ylab="Variance in Coverage")

#assembled together
wilcox.test(bo.cov$CovVariance,bd.cov$CovVariance,paired=T,"less")


###Heterozygosity
#assembled separately
hist(o.cov$PropHet,col=alpha("cadetblue",0.5))
hist(d.cov$PropHet,col=alpha("coral1",0.5),add=T)
boxplot(o.cov$PropHet,d.cov$PropHet,col=c("cadetblue","coral1"),
        names=c("sdRAD","ddRAD"),ylab="Proportion Heterozygous")
text(x=c(1,2),y=c(21,21),c("250425","69109"))
wilcox.test(o.cov$PropHet,d.cov$PropHet,"less")

#what about in 30 males 30 females from ddRAD and the sdRAD?
d.mal<-colnames(drad)[grep("PRM",colnames(drad))]
d.fem<-colnames(drad)[grep("FEM",colnames(drad))]
dm.sub<-sample(d.mal,30)
df.sub<-sample(d.fem,30)
d.sub<-cbind(drad[,1:9],drad[,dm.sub],drad[,df.sub])
ds.cov<-do.call("rbind",apply(d.sub,1,vcf.cov.loc,subset=c(dm.sub,df.sub)))

wilcox.test(o.cov$PropHet,ds.cov$PropHet,"less")
wilcox.test(o.cov$AvgCovRatio,ds.cov$AvgCovRatio,"less")
wilcox.test(o.cov$CovVariance,ds.cov$CovVariance,"greater")

#assembled together
wilcox.test(bo.cov$PropHet,bd.cov$PropHet,paired=T,"less")


lc.comp<-data.frame(LibraryPrep=c(rep("sdRAD",nrow(o.cov)),rep("ddRAD",nrow(d.cov)),
                                  rep("sdRAD",nrow(bo.cov)),rep("ddRAD",nrow(bd.cov))), 
                    Assembly=c(rep("Alone",nrow(o.cov)),rep("Alone",nrow(d.cov)),
                               rep("Together",nrow(bo.cov)),rep("Together",nrow(bd.cov))),
                    CovVariance=c(o.cov$CovVariance,d.cov$CovVariance,bo.cov$CovVariance,bd.cov$CovVariance),
                    PropHet=c(o.cov$PropHet,d.cov$PropHet,bo.cov$PropHet,bd.cov$PropHet))

summary(aov(log(lc.comp$CovVariance+1)~lc.comp$LibraryPrep*lc.comp$Assembly))
summary(aov(lc.comp$PropHet~lc.comp$LibraryPrep*lc.comp$Assembly))

jpeg("VarianceInCov.jpeg",height=8.5,width=7,units="in",res=300)
par(mfrow=c(2,1),oma=c(1,2,1,1),mar=c(2,2,2,2))
boxplot(log(lc.comp$CovVariance+1)~lc.comp$LibraryPrep*lc.comp$Assembly,col=c("cadetblue","coral1"),
        names=F,xaxt='n')
axis(1,at=c(1.5,3.5),c("Alone","Together"))
mtext("log(Variance in Coverage)",2,outer=F,line=2)
legend("topleft",ncol=2,c("sdRAD","ddRAD"),pch=15,col=c("cadetblue","coral1"),bty='n')
boxplot(lc.comp$PropHet~lc.comp$LibraryPrep*lc.comp$Assembly,col=c("cadetblue","coral1"),
        names=F,xaxt='n')
mtext("Proportion Heterozygotes",2,outer=F,line=2)
axis(1,at=c(1.5,3.5),c("Alone","Together"))
dev.off()

###########################LOOKING FOR THE SAME LOCI#################################
#drad,orad,both
drad$SNP<-paste(drad$`#CHROM`,drad$POS,sep=".")
orad$SNP<-paste(orad$`#CHROM`,orad$POS,sep=".")
both$SNP<-paste(both$`#CHROM`,both$POS,sep=".")
locus.info<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
d.cov$SNP<-paste(d.cov$Chrom,d.cov$Pos,sep=".")
o.cov$SNP<-paste(o.cov$Chrom,o.cov$Pos,sep=".")
b.cov$SNP<-paste(b.cov$Chrom,b.cov$Pos,sep=".")
#sdRAD-ddRAD
od.loci<-merge(orad[,locus.info],drad[,locus.info],"SNP")
dim(od.loci[(od.loci$REF.x != od.loci$REF.y & od.loci$REF.x != od.loci$ALT.y) |
               (od.loci$REF.x != od.loci$REF.y & od.loci$REF.y != od.loci$ALT.x),])
    
               dim(od.loci[(od.loci$REF.x != od.loci$REF.y & od.loci$ALT.x != od.loci$ALT.y &
                  od.loci$REF.x != od.loci$ALT.y & od.loci$REF.y != od.loci$ALT.x) ,])
od.cov<-merge(o.cov,d.cov,"SNP")

wilcox.test(od.cov$PropHet.x,od.cov$PropHet.y,paired=T,"less")
wilcox.test(od.cov$CovVariance.x,od.cov$CovVariance.y,paired=T,"greater")
wilcox.test(od.cov$AvgCovRatio.x,od.cov$AvgCovRatio.y,paired=T,"greater")

#sdRAD-both
ob.loci<-merge(orad[,locus.info],both[,locus.info],"SNP")
dim(ob.loci[(ob.loci$REF.x != ob.loci$REF.y & ob.loci$REF.x != ob.loci$ALT.y) |
              (ob.loci$REF.x != ob.loci$REF.y & ob.loci$REF.y != ob.loci$ALT.x),])
dim(ob.loci[(ob.loci$REF.x != ob.loci$REF.y & ob.loci$ALT.x != ob.loci$ALT.y &
               ob.loci$REF.x != ob.loci$ALT.y & ob.loci$REF.y != ob.loci$ALT.x) ,])
ob.cov<-merge(o.cov,b.cov,"SNP")
wilcox.test(ob.cov$PropHet.x,ob.cov$PropHet.y,paired=T,"less")
wilcox.test(ob.cov$CovVariance.x,ob.cov$CovVariance.y,paired=T,"greater")
wilcox.test(ob.cov$AvgCovRatio.x,ob.cov$AvgCovRatio.y,paired=T,"less")

#ddRAD-both
db.loci<-merge(drad[,locus.info],both[,locus.info],"SNP")
dim(db.loci[(db.loci$REF.x != db.loci$REF.y & db.loci$REF.x != db.loci$ALT.y) |
              (db.loci$REF.x != db.loci$REF.y & db.loci$REF.y != db.loci$ALT.x),])
dim(db.loci[(db.loci$REF.x != db.loci$REF.y & db.loci$ALT.x != db.loci$ALT.y &
               db.loci$REF.x != db.loci$ALT.y & db.loci$REF.y != db.loci$ALT.x) ,])
db.cov<-merge(d.cov,b.cov,"SNP")
wilcox.test(db.cov$PropHet.x,db.cov$PropHet.y,paired=T,"greater")
wilcox.test(db.cov$CovVariance.x,db.cov$CovVariance.y,paired=T,"less")
wilcox.test(db.cov$AvgCovRatio.x,db.cov$AvgCovRatio.y,paired=T,"less")

#sd/dd-both
sdb.loci<-merge(od.loci,both[,locus.info],"SNP")
dim(sdb.loci[(sdb.loci$REF.x != sdb.loci$REF.y & sdb.loci$REF.x != sdb.loci$REF
              & sdb.loci$REF.x != sdb.loci$ALT.y & sdb.loci$REF.x != sdb.loci$ALT) |
               (sdb.loci$ALT.x != sdb.loci$REF.y & sdb.loci$ALT.x != sdb.loci$REF
                & sdb.loci$ALT.x != sdb.loci$ALT.y & sdb.loci$ALT.x != sdb.loci$ALT) |
               (sdb.loci$REF.y != sdb.loci$REF.x & sdb.loci$REF.y != sdb.loci$REF
                & sdb.loci$REF.y != sdb.loci$ALT.x & sdb.loci$REF.y != sdb.loci$ALT) |
               (sdb.loci$ALT.y != sdb.loci$REF.x & sdb.loci$ALT.y != sdb.loci$REF
                & sdb.loci$ALT.y != sdb.loci$ALT.y & sdb.loci$ALT.y != sdb.loci$ALT) |
               (sdb.loci$REF != sdb.loci$REF.x & sdb.loci$REF != sdb.loci$REF.x
                & sdb.loci$REF != sdb.loci$ALT.y & sdb.loci$REF != sdb.loci$ALT.x) |
               (sdb.loci$ALT != sdb.loci$REF.y & sdb.loci$ALT != sdb.loci$REF.x
                & sdb.loci$ALT != sdb.loci$ALT.y & sdb.loci$ALT.x != sdb.loci$ALT) ,])

#########################FSTS FROM INDEPENDENT ANALYSES##############################
pop1<-colnames(orad)[10:ncol(orad)]
pop2<-colnames(drad)[10:ncol(drad)]
o.share<-orad[orad$SNP %in% od.loci$SNP,]
d.share<-drad[drad$SNP %in% od.loci$SNP,]
od.vcf<-merge(orad,drad,"SNP")


od.fst<-do.call("rbind",apply(o.share,1,fst.two.vcf,vcf2=d.share,match.index="SNP"))
write.table(od.fst,"orad-drad.fst.txt",col.names=T,row.names=F,quote=F,sep='\t')

#FROM SAME ANALYSIS
both$SNP<-paste(both$`#CHROM`,both$POS,sep=".")
both.o<-cbind(both[,locus.info],both[,o.ind])
both.d<-cbind(both[,locus.info],both[,d.ind])
od.both.fst<-do.call("rbind",apply(both.o,1,fst.two.vcf,vcf2=both.d,match.index="SNP"))
od.both.fst$SNP<-paste(od.both.fst$Chrom,od.both.fst$Pos,sep=".")
od.fst.1<-od.both.fst[od.both.fst$Fst ==1,]
bd.cov$SNP<-paste(bd.cov$Chrom,bd.cov$Pos,sep=".")
bo.cov$SNP<-paste(bo.cov$Chrom,bo.cov$Pos,sep=".")
bd.cov.pass<-bd.cov[bd.cov$AvgCovTotal > 3 & bd.cov$AvgCovTotal <=50,"SNP"]
bo.cov.pass<-bo.cov[bo.cov$AvgCovTotal > 3 & bo.cov$AvgCovTotal <=50,"SNP"]
cov.pass<-bd.cov.pass[bd.cov.pass %in% bo.cov.pass]
plot(od.both.fst[od.both.fst$SNP %in% cov.pass,"Fst"])

#PLOTTING
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
scaffs<-levels(as.factor(both[,1]))
scaffs[1:22]<-lgs

jpeg("sd-dd_Fst.jpeg",height=10,width=7.5,units="in",res=300)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(1,2,1,1))
od<-fst.plot(od.fst[!is.na(od.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
  fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
  groups=as.factor(scaffs[scaffs %in% levels(factor(od.fst$Chrom))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(od[od$Chrom ==lgs[i],"Pos"]),y=-1,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(od[od$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately",3,cex=0.75)
odb<-fst.plot(od.both.fst[od.both.fst$Num1>30 & od.both.fst$Num2 > 192,],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
              groups=as.factor(scaffs[scaffs %in% levels(factor(od.both.fst[od.both.fst$Num1>30 & od.both.fst$Num2 > 192,"Chrom"]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odb[odb$Chrom ==lgs[i],"Pos"]),y=-1.7,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odb[odb$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Together",3,cex=0.75)
odbs<-fst.plot(od.both.fst[od.both.fst$SNP %in% od.loci$SNP & od.both.fst$Num1>30 & od.both.fst$Num2 > 192,],
  ci.dat=c(-10,10),sig.col=c("black","black"), 
  fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
  groups=as.factor(scaffs[scaffs %in% 
    levels(factor(od.both.fst[od.both.fst$SNP %in% od.loci$SNP & od.both.fst$Num1>30 & od.both.fst$Num2 > 192,"Chrom"]))]))
last<-0
for(i in 1:length(lgs)){
  if(lgs[i] %in% levels(factor(odbs$Chrom))){
  text(x=mean(odbs[odbs$Chrom ==lgs[i],"Pos"]),y=-0.66,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odbs[odbs$Chrom ==lgs[i],"Pos"])
  }
}
mtext("sdRAD-ddRAD Analyzed Together, Only Loci In Separate Analyses",3,cex=0.75)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)

odbc<-fst.plot(od.both.fst[od.both.fst$SNP %in% cov.pass & !is.na(od.both.fst$Fst),],
               ci.dat=c(-10,10),sig.col=c("black","black"), 
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                levels(factor(od.both.fst[od.both.fst$SNP %in% cov.pass & !is.na(od.both.fst$Fst),"Chrom"]))]))
last<-0
for(i in 1:length(lgs)){
  if(lgs[i] %in% levels(factor(odbc$Chrom))){
    text(x=mean(odbc[odbc$Chrom ==lgs[i],"Pos"]),y=-1.5,
         labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
    last<-max(odbc[odbc$Chrom ==lgs[i],"Pos"])
  }
}
mtext("sdRAD-ddRAD Analyzed Together, High Coverage Removed",3,cex=0.75)

dev.off()

fst.aov.dat<-data.frame(Fst=c(od$Fst,odb$Fst,odbs$Fst),
  Type=c(rep("Alone",nrow(od)),rep("Together",nrow(odb)),rep("Together,Alone Loci",nrow(odbs))))
fst.aov<-aov(Fst~Type,dat=fst.aov.dat)

##############################DRAD DIFFERENT PLATES##################################

drad$SNP<-paste(drad$`#CHROM`,drad$ID,drad$POS,sep=".")

d.cov<-do.call("rbind",apply(drad,1,vcf.cov.loc,subset=d.ind))
d.cov$SNP<-paste(d.cov$Chrom,d.cov$Locus,d.cov$Pos,sep=".")
d.cov$AvgCovRatio[d.cov$AvgCovRatio=="Inf"]<-0
d.keep<-d.cov[d.cov$PropMissing>0.5,]
plate1<-unlist(as.list(read.delim("../plate1.txt",stringsAsFactors = F)))
plate2<-unlist(as.list(read.delim("../plate2.txt",stringsAsFactors = F)))
plate3<-unlist(as.list(read.delim("../plate3.txt",stringsAsFactors = F)))
plate4<-unlist(as.list(read.delim("../plate4.txt",stringsAsFactors = F)))

d1<-cbind(drad[drad$SNP %in% d.keep$SNP,1:9],drad$SNP[drad$SNP %in% d.keep$SNP],drad[drad$SNP %in% d.keep$SNP,colnames(drad) %in% plate1])
d2<-cbind(drad[drad$SNP %in% d.keep$SNP,1:9],drad$SNP[drad$SNP %in% d.keep$SNP],drad[drad$SNP %in% d.keep$SNP,colnames(drad) %in% plate2])
d3<-cbind(drad[drad$SNP %in% d.keep$SNP,1:9],drad$SNP[drad$SNP %in% d.keep$SNP],drad[drad$SNP %in% d.keep$SNP,colnames(drad) %in% plate3])
d4<-cbind(drad[drad$SNP %in% d.keep$SNP,1:9],drad$SNP[drad$SNP %in% d.keep$SNP],drad[drad$SNP %in% d.keep$SNP,colnames(drad) %in% plate4])
p1<-grep("sample",colnames(d1))
p2<-grep("sample",colnames(d2))
p3<-grep("sample",colnames(d3))
p4<-grep("sample",colnames(d4))
d1.cov<-do.call("rbind",apply(d1,1,vcf.cov.loc,subset=p1))
d2.cov<-do.call("rbind",apply(d2,1,vcf.cov.loc,subset=p2))
d3.cov<-do.call("rbind",apply(d3,1,vcf.cov.loc,subset=p3))
d4.cov<-do.call("rbind",apply(d4,1,vcf.cov.loc,subset=p4))
d1.cov$Plate<-"plate1"
d2.cov$Plate<-"plate2"
d3.cov$Plate<-"plate3"
d4.cov$Plate<-"plate4"
dcomp<-rbind(d1.cov,d2.cov,d3.cov,d4.cov)
dcomp$AvgCovRatio[dcomp$AvgCovRatio=="Inf"]<-0
summary(aov(dcomp$PropHet~dcomp$Plate))
TukeyHSD(aov(dcomp$PropHet~dcomp$Plate))

summary(aov(dcomp$AvgCovTotal~dcomp$Plate))
TukeyHSD(aov(dcomp$AvgCovTotal~dcomp$Plate))

summary(aov(dcomp$AvgCovRatio~dcomp$Plate))
TukeyHSD(aov(dcomp$AvgCovRatio~dcomp$Plate))
drl<-lm(dcomp$AvgCovRatio~dcomp$Locus*dcomp$Plate)
plot(dcomp$AvgCovRatio~dcomp$Locus,col=alpha(as.numeric(as.factor(dcomp$Plate)),0.5),pch=19)
legend("top",ncol=2,levels(factor(dcomp$Plate)),col=alpha(as.numeric(as.factor(levels(as.factor(dcomp$Plate)))),0.5),pch=19,bty='n')
summary(aov(dcomp$CovVariance~dcomp$Plate))
TukeyHSD(aov(dcomp$CovVariance~dcomp$Plate))
dvl<-lm(dcomp$CovVariance~dcomp$Locus*dcomp$Plate)
plot(dcomp$CovVariance~dcomp$Locus,col=alpha(as.numeric(as.factor(dcomp$Plate)),0.5),pch=19)
legend("top",ncol=2,levels(factor(dcomp$Plate)),col=alpha(as.numeric(as.factor(levels(as.factor(dcomp$Plate)))),0.5),pch=19,bty='n')

summary(aov(dcomp$PropMissing~dcomp$Plate))
TukeyHSD(aov(dcomp$PropMissing~dcomp$Plate))

############################ORIGINAL INVESTIGATION###################################
setwd("E:/ubuntushare/SCA/results/biallelic/both_datasets/")
summary<-read.delim("gwsca_summary_datasets.txt")

ddRAD<-summary[summary$Pop=="ddRAD",]
oRAD<-summary[summary$Pop=="oRAD",]
colB<-c(rep("dd",nrow(ddRAD)),rep("o",nrow(oRAD)))
colA<-c(ddRAD$Allele1Freq,oRAD$Allele1Freq)
dt<-data.frame(colA,colB)

png("AlleleFreq_twotypes.png",height=7,width=7,units="in",res=300)
ggplot(dt, aes(x = colA, fill = colB)) + 
  geom_histogram(col = "black", alpha = 0.5, 
                 position = "identity") +
  scale_fill_discrete(name = "Library Type",
                      labels = c("ddRAD","sdRAD"))
dev.off()

fst<-read.delim("gwsca_fsts_datasets.txt")
mean(fst$ddRAD.oRAD[fst$ddRAD.oRAD>=0]) #0.0616421
png("DDvsSD_allFsts.png",height=7,width=7,units="in",res=300)
plot(fst$ddRAD.oRAD[fst$ddRAD.oRAD>=0], xaxt='n',ylab="Fst",pch=19) #0.0616421
dev.off()


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


