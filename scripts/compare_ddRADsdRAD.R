#Author: Sarah P. Flanagan
#Last updated: 13 September 2016
#Date Started: 24 February 2016
#Purpose: compare RAD libraries prepared by two different methods

rm(list=ls())
library(ggplot2)
setwd("~/Projects/SCA/results")

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
    }))))/pres
  het<-unlist(lapply(vcf.row[subset],function(x){ 
    strsplit(as.character(x),split=":")[[1]][1]
  }))
  het<-length(het[het=="0/1" | het=="1/0"])
  return(data.frame(Locus=vcf.row["ID"],NumMissing=miss, NumPresent=pres,AvgCovRef=ref,
    AvgCovAlt=alt, AvgCovTotal=tot, NumHet=het,stringsAsFactors = F))
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
  
######FILES#####
drad<-parse.vcf("drad.vcf")
orad<-parse.vcf("orad.vcf")
both<-parse.vcf("both.vcf")

#########ANALYSIS##########
#coverage per locus
o.ind<-grep("orad",colnames(both),value=T)
d.ind<-grep("sample",colnames(both),value=T)

o.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=o.ind))
d.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=d.ind))
loc.cov<-merge(o.cov,d.cov,by="Locus")
loc.cov$oRAD.AvgCovRatio<-loc.cov$AvgCovRef.x/loc.cov$AvgCovAlt.x
loc.cov$dRAD.AvgCovRatio<-loc.cov$AvgCovRef.y/loc.cov$AvgCovAlt.y
#because some of them have no alt
loc.cov$oRAD.AvgCovRatio[loc.cov$oRAD.AvgCovRatio=="Inf"]<-loc.cov[loc.cov$oRAD.AvgCovRatio=="Inf","AvgCovRef.x"]
loc.cov$dRAD.AvgCovRatio[loc.cov$dRAD.AvgCovRatio=="Inf"]<-loc.cov[loc.cov$dRAD.AvgCovRatio=="Inf","AvgCovRef.y"]
loc.cov$oRAD.PropMiss<-loc.cov$NumMissing.x/(loc.cov$NumMissing.x+loc.cov$NumPresent.x)
loc.cov$dRAD.PropMiss<-loc.cov$NumMissing.y/(loc.cov$NumMissing.y+loc.cov$NumPresent.y)
loc.cov$oRAD.PropHet<-loc.cov$NumHet.x/loc.cov$NumPresent.x
loc.cov$dRAD.PropHet<-loc.cov$NumHet.y/loc.cov$NumPresent.y
t.test(loc.cov$oRAD.AvgCovRatio,loc.cov$dRAD.AvgCovRatio,paired=T,alternative="less")
#private alleles
dim(loc.cov[(loc.cov$AvgCovRef.x == 0 & loc.cov$AvgCovAlt.y == 0) | (loc.cov$AvgCovAlt.x == 0 & loc.cov$AvgCovRef.y == 0),])#2042
dim(loc.cov[loc.cov$NumHet.x==0,])#72307
dim(loc.cov[loc.cov$NumHet.y==0,])#14436
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


