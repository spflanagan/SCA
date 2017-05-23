#Author: Sarah P. Flanagan
#Last updated: 23 May 2017
#Date Started: 24 February 2016
#Purpose: compare RAD libraries prepared by two different methods

rm(list=ls())
library(ggplot2)
setwd("~/Projects/SCA/results")
source("../../gwscaR/R/gwscaR.R")
#################FUNCTIONS####################
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
  
fst.two.vcf<-function(vcf1.row,vcf2,match.index, cov.thresh=0.2){
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
      num.ind<-length(gt1)-10
      gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
      gt1[gt1=="1/0"]<-"0/1"
      gt1<-gsub(pattern = "0",replacement = vcf1.row["REF"],gt1)
      gt1<-gsub(pattern = "1",replacement = vcf1.row["ALT"],gt1)
      if(length(gt1)/num.ind>=cov.thresh){
        al1<-unlist(strsplit(as.character(gt1),split = "/"))
        gt2<-unlist(lapply(vcf2.row,function(x){ 
          c<-strsplit(as.character(x),split=":")[[1]][1]
          return(c)
        }))
        num.ind<-length(gt2)
        gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
        gt2[gt2=="1/0"]<-"0/1"
        gt2<-gsub(pattern = "0",replacement = vcf2.row["REF"],gt2)
        gt2<-gsub(pattern = "1",replacement = vcf2.row["ALT"],gt2)
        if(length(gt2)/num.ind>=cov.thresh){
          al2<-unlist(strsplit(as.character(gt2),split="/"))
           #calculate frequencies
          freq1<-summary(factor(al1))/sum(summary(factor(al1)))	
          freq2<-summary(factor(al2))/sum(summary(factor(al2)))	
          freqall<-summary(as.factor(c(al1,al2)))/
            sum(summary(as.factor(c(al1,al2))))
          hets<-c(names(freq1)[2],names(freq2)[2])
          if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
            hs1<-2*freq1[1]*freq1[2]
            hs2<-2*freq2[1]*freq2[2]
            hs<-mean(c(hs1,hs2))
            ht<-2*freqall[1]*freqall[2]
            fst<-(ht-hs)/ht
          } else {
            hs1<-1-sum(freq1*freq1)
            hs2<-1-sum(freq2*freq2)
            if(length(freqall)<=1){ fst<-0 }
            else{ 
              ht<-2*freqall[1]*freqall[2]
              fst<-NA
            }
          }
        }
        else {
          fst<-NA #gt2 doesn't pass coverage threshold
        }
      }else {
        fst<-NA #it doesn't pass the coverage threshold
      }
    }#end else good to go
  }#end else vcf2

  return(data.frame(Chrom=vcf1.row["#CHROM"],Pos=vcf1.row["POS"],
    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
    Num1=length(gt1),Num2=(length(gt2))))
}#end function


calc.afs.vcf<-function(vcf.row){
  #use in conjunction with apply
  #e.g. apply(vcf,1,afs.vcf)
  gt1<-unlist(lapply(vcf.row,function(x){ 
    c<-strsplit(as.character(x),split=":")[[1]][1]
    return(c)
  }))
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1=="1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  al1<-unlist(strsplit(as.character(gt1),split = "/"))
  #calculate frequencies
  freq1<-summary(factor(al1))/sum(summary(factor(al1)))	
  if(length(freq1)==1)
  {
    if(names(freq1)==vcf.row["REF"])
    {
      freq1<-c(freq1,0)
      names(freq1)<-unlist(c(vcf.row["REF"],vcf.row["ALT"]))
    }
    else
    {
      freq1<-c(freq1,0)
      names(freq1)<-unlist(c(vcf.row["ALT"],vcf.row["REF"]))
    }
  }
  return(data.frame(Chrom=vcf.row["#CHROM"], Pos=vcf.row["POS"], Ref=vcf.row["REF"],
    RefFreq=freq1[names(freq1) %in% vcf.row["REF"]],
    Alt=vcf.row["ALT"],AltFreq=freq1[names(freq1) %in% vcf.row["ALT"]]))
}

fst.one.vcf<-function(vcf.row,group1,group2, cov.thresh=0.2){
  
  hs1<-hs2<-hs<-ht<-0
  freqall<-gt1<-gt2<-NULL
  gt1<-unlist(lapply(vcf.row[group1],function(x){ 
    c<-strsplit(as.character(x),split=":")[[1]][1]
    return(c)
  }))
  num.ind<-length(gt1)
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1=="1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  if(length(gt1)/num.ind>=cov.thresh){
    al1<-unlist(strsplit(as.character(gt1),split = "/"))
    gt2<-unlist(lapply(vcf.row[group2],function(x){ 
      c<-strsplit(as.character(x),split=":")[[1]][1]
      return(c)
    }))
    num.ind<-length(gt2)
    gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
    gt2[gt2=="1/0"]<-"0/1"
    gt2<-gsub(pattern = "0",replacement = vcf.row["REF"],gt2)
    gt2<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt2)
    if(length(gt2)/num.ind>=cov.thresh){
      al2<-unlist(strsplit(as.character(gt2),split="/"))
      #calculate frequencies
      freq1<-summary(factor(al1))/sum(summary(factor(al1)))	
      freq2<-summary(factor(al2))/sum(summary(factor(al2)))	
      freqall<-summary(as.factor(c(al1,al2)))/
        sum(summary(as.factor(c(al1,al2))))
      hets<-c(names(freq1)[2],names(freq2)[2])
      if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
        hs1<-2*freq1[1]*freq1[2]
        hs2<-2*freq2[1]*freq2[2]
        hs<-mean(c(hs1,hs2))
        ht<-2*freqall[1]*freqall[2]
        fst<-(ht-hs)/ht
      } else {
        hs1<-1-sum(freq1*freq1)
        hs2<-1-sum(freq2*freq2)
        if(length(freqall)<=1){ fst<-0 }
        else{ 
          ht<-2*freqall[1]*freqall[2]
          fst<-NA
        }
      }
    }
    else {
      fst<-NA #gt2 doesn't pass coverage threshold
    }
  }else {
    fst<-NA #it doesn't pass the coverage threshold
  }

  return(data.frame(Chrom=vcf.row["#CHROM"],Pos=vcf.row["POS"],
                  Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
                  Num1=length(gt1),Num2=length(gt2)))
}#end function fst.one.vcf

choose.one.snp<-function(vcf){
  keep.col<-colnames(vcf)
  vcf$id.pos<-paste(vcf$ID,vcf$POS,sep=".")
  sub.vcf<-tapply(vcf$id.pos,vcf$ID, sample,size=1)
  new.vcf<-vcf[vcf$id.pos %in% sub.vcf,keep.col]
  return(new.vcf)
}



fst.one.plink<-function(raw,group1, group2, cov.thresh=0.2){
  fst.dat<-data.frame(Locus=character(),
             Hs1=numeric(),Hs2=numeric(),Hs=numeric(),Ht=numeric(),Fst=numeric(),NumAlleles=numeric(),
             Num1=numeric(),Num2=numeric(),stringsAsFactors=F)
  grp1<-raw[raw$IID %in% group1,]
  grp2<-raw[raw$IID %in% group2,]
  for(i in 7:ncol(raw)){
    na1<-length(grp1[is.na(grp1[,i]),i])/nrow(grp1)
    na2<-length(grp2[is.na(grp2[,i]),i])/nrow(grp2)
    gt1<-grp1[!is.na(grp1[,i]),i]
    gt2<-grp2[!is.na(grp2[,i]),i]
    gt1[gt1=="1"]<-"1/2"
    gt1[gt1=="2"]<-"2/2"
    gt1[gt1=="0"]<-"1/1"
    gt2[gt2=="1"]<-"1/2"
    gt2[gt2=="2"]<-"2/2"
    gt2[gt2=="0"]<-"1/1"
    
    if(na1<=(1-cov.thresh)){
      al1<-unlist(strsplit(as.character(gt1),split = "/"))
      if(na2<=(1-cov.thresh)){
        al2<-unlist(strsplit(as.character(gt2),split="/"))
        #calculate frequencies
        freq1<-summary(factor(al1))/sum(summary(factor(al1)))	
        freq2<-summary(factor(al2))/sum(summary(factor(al2)))	
        freqall<-summary(as.factor(c(al1,al2)))/
          sum(summary(as.factor(c(al1,al2))))
        if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
          hs1<-2*freq1[1]*freq1[2]
          hs2<-2*freq2[1]*freq2[2]
          hs<-mean(c(hs1,hs2))
          ht<-2*freqall[1]*freqall[2]
          fst<-(ht-hs)/ht
        } else {
          hs1<-1-sum(freq1*freq1)
          hs2<-1-sum(freq2*freq2)
          if(length(freqall)<=1){ fst<-0 }
          else{ 
            ht<-2*freqall[1]*freqall[2]
            fst<-NA
          }
        }
      }
      else {
        fst<-NA #gt2 doesn't pass coverage threshold
      }
    }else {
      fst<-NA #it doesn't pass the coverage threshold
    }
    fst.dat[(i-6),]<-cbind(as.character(colnames(raw)[i]),hs1,hs2,as.numeric(hs),ht,fst,length(freqall),length(gt1),length(gt2))
  }
  return(fst.dat)
}#end fst.one.plink

#Author: Sarah P. Flanagan
#Last updated: 23 September 2016
#Date Started: 24 February 2016
#Purpose: compare RAD libraries prepared by two different methods

rm(list=ls())
library(ggplot2)
setwd("~/Projects/SCA/results")
source("../../gwscaR/R/gwscaR.R")

#################SET COLORS#################
ddtog.col<-"#f4a582"
ddsep.col<-"#ca0020"
sdtog.col<-"#92c5de"
sdsep.col<-"#0571b0"
ddsdtog.col<-"#807dba" #or "#6a51a3"
ddsdsep.col<-"#3f007d"
sexsel.col<-"navyblue"
viasel.col<-"goldenrod"
#################FILES#################
#Original files
#drad<-parse.vcf("drad.vcf")
#orad<-parse.vcf("orad.vcf")
#both<-parse.vcf("both.vcf")
#get snp subsets and write to file
#both.sub<-choose.one.snp(both)
#write.table(both.sub,"both.sub.vcf",quote=F,sep="\t",col.names = T,row.names=F)
#drad.sub<-choose.one.snp(drad)
#write.table(drad.sub,"drad.sub.vcf",quote=F,sep="\t",col.names = T,row.names=F)
#orad.sub<-choose.one.snp(orad)
#write.table(orad.sub,"orad.sub.vcf",quote=F,sep="\t",col.names = T,row.names=F)
drad<-parse.vcf("drad.sub.vcf")
orad<-parse.vcf("orad.sub.vcf")
both<-parse.vcf("both.sub.vcf")
drad$SNP<-paste(drad$`#CHROM`,as.numeric(drad$POS),sep=".")
orad$SNP<-paste(orad$`#CHROM`,as.numeric(orad$POS),sep=".")
both$SNP<-paste(both$`#CHROM`,as.numeric(both$POS),sep=".")
locus.info<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
o.ind<-grep("orad",colnames(both),value=T)
d.ind<-grep("sample",colnames(both),value=T)



#plink subset files - don't use
#dsub<-read.delim("drad.subset.raw",sep=" ")
#dmap<-read.delim("drad.map",header=F)
#osub<-read.delim("orad.subset.raw",sep=" ")
#omap<-read.delim("orad.map",header=F)
#bsub<-read.delim("both.subset.raw",sep=" ")
#bmap<-read.delim("both.map",header=F)

#dmap$SNP<-paste(dmap$V1,dmap$V4,sep=".")
#omap$SNP<-paste(omap$V1,omap$V4,sep=".")
#bmap$SNP<-paste(bmap$V1,bmap$V4,sep=".")
#colnames(omap)<-c("Chrom","oLocus","dist","Pos","SNP")
#colnames(dmap)<-c("Chrom","dLocus","dist","Pos","SNP")
#colnames(bmap)<-c("Chrom","Locus","dist","Pos","SNP")
#################ANALYSIS#################
##ran these already--read in files##
#bo.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=o.ind))
#bd.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=d.ind))
#b.cov<-do.call("rbind",apply(both,1,vcf.cov.loc,subset=c(o.ind,d.ind)))
#bo.cov$SNP<-paste(bo.cov$Chrom,bo.cov$Pos,sep=".")
#bd.cov$SNP<-paste(bd.cov$Chrom,bd.cov$Pos,sep=".")
#b.cov$SNP<-paste(b.cov$Chrom,b.cov$Pos,sep=".")
#o.cov<-do.call("rbind",apply(orad,1,vcf.cov.loc,subset=o.ind))
#o.cov$SNP<-paste(o.cov$Chrom,o.cov$Pos,sep=".")
#d.cov<-do.call("rbind",apply(drad,1,vcf.cov.loc,subset=d.ind))
#d.cov$SNP<-paste(d.cov$Chrom,d.cov$Pos,sep=".")
#loc.cov<-merge(o.cov,d.cov,by="SNP")
#because some of them have no alt
#loc.cov$AvgCovRatio.x[loc.cov$AvgCovRatio.x=="Inf"]<-loc.cov[loc.cov$AvgCovRatio.x=="Inf","AvgCovRef.x"]
#loc.cov$AvgCovRatio.y[loc.cov$AvgCovRatio.y=="Inf"]<-loc.cov[loc.cov$AvgCovRatio.y=="Inf","AvgCovRef.y"]
#write.csv(b.cov,"Both_overallCov.csv")
#write.csv(bo.cov,"Both_sdRADCov.csv")
#write.csv(bd.cov,"Both_ddRADCov.csv")
#write.csv(o.cov,"Sep_sdRADCov.csv")
#write.csv(d.cov,"Sep_ddRADCov.csv")
#write.csv(loc.cov,"Sep_Combined.csv")
b.cov<-read.csv("Both_overallCov.csv",row.names=1)
b.cov$SNP<-paste(b.cov$Chrom,b.cov$Pos,sep=".")
bo.cov<-read.csv("Both_sdRADCov.csv",row.names=1)
bo.cov$SNP<-paste(bo.cov$Chrom,bo.cov$Pos,sep=".")
bd.cov<-read.csv("Both_ddRADCov.csv",row.names=1)
bd.cov$SNP<-paste(bd.cov$Chrom,bd.cov$Pos,sep=".")
o.cov<-read.csv("Sep_sdRADCov.csv",row.names=1)
o.cov$SNP<-paste(o.cov$Chrom,o.cov$Pos,sep=".")
d.cov<-read.csv("Sep_ddRADCov.csv",row.names=1)
d.cov$SNP<-paste(d.cov$Chrom,d.cov$Pos,sep=".")
loc.cov<-read.csv("Sep_Combined.csv",row.names=1)
loc.cov$SNP<-paste(loc.cov$Chrom.x,loc.cov$Pos.x,sep=".")

bo.cov$AvgCovRatio[is.na(bo.cov$AvgCovRatio)]<-0
bo.cov$AvgCovRatio[bo.cov$AvgCovRatio=="Inf"]<-bo.cov[bo.cov$AvgCovRatio=="Inf","AvgCovRef"]
bd.cov$AvgCovRatio[bd.cov$AvgCovRatio=="Inf"]<-bd.cov[bd.cov$AvgCovRatio=="Inf","AvgCovRef"]

##SUBSET
b.cov<-b.cov[b.cov$SNP %in% both$SNP,]
bo.cov<-bo.cov[bo.cov$SNP %in% both$SNP & bo.cov$SNP %in% orad$SNP,]
bd.cov<-bd.cov[bd.cov$SNP %in% both$SNP & bd.cov$SNP %in% drad$SNP,]
o.cov<-o.cov[o.cov$SNP %in% orad$SNP,]
d.cov<-d.cov[d.cov$SNP %in% drad$SNP,]
loc.cov<-loc.cov[loc.cov$SNP %in% orad$SNP & loc.cov$SNP %in% drad$SNP,]
####Average Coverage####
#compare sd vs dd
lc.comp<-data.frame(LibraryPrep=c(rep("sdRAD",nrow(o.cov)),rep("ddRAD",nrow(d.cov)),
    rep("sdRAD",nrow(bo.cov)),rep("ddRAD",nrow(bd.cov))), 
  Assembly=c(rep("Alone",nrow(o.cov)),rep("Alone",nrow(d.cov)),
             rep("Together",nrow(bo.cov)),rep("Together",nrow(bd.cov))),
  AvgCovRatio=c(o.cov$AvgCovRatio,d.cov$AvgCovRatio,bo.cov$AvgCovRatio,bd.cov$AvgCovRatio),
  AvgCovTotal=c(o.cov$AvgCovTotal,d.cov$AvgCovTotal,bo.cov$AvgCovTotal,bd.cov$AvgCovTotal))
#NumReadsPerInd=c(oc$NumReadsAlone,dc$NumReadsAlone,oc$NumReadsTogether,dc$NumReadsTogether)

#total coverage
summary(aov(log(lc.comp$AvgCovTotal)~lc.comp$LibraryPrep*lc.comp$Assembly))
TukeyHSD(aov(log(lc.comp$AvgCovTotal)~lc.comp$LibraryPrep*lc.comp$Assembly))#use this!
summary(lm(log(lc.comp$AvgCovTotal)~lc.comp$LibraryPrep*lc.comp$Assembly))
boxplot(log(lc.comp$AvgCovTotal)~lc.comp$LibraryPrep*lc.comp$Assembly)
#ref/alt ratio
summary(aov(lc.comp$AvgCovRatio~lc.comp$LibraryPrep*lc.comp$Assembly))
TukeyHSD(aov(log(lc.comp$AvgCovRatio+1)~lc.comp$LibraryPrep*lc.comp$Assembly))#use this!
summary(lm(lc.comp$AvgCovRatio~lc.comp$LibraryPrep*lc.comp$Assembly))
boxplot(log(lc.comp$AvgCovRatio+1)~lc.comp$LibraryPrep*lc.comp$Assembly)

####PLOT: Coverage comparison boxplot####
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

####Coverage by individual--assembled in one####
##Generate the dataset##
#o.icov<-as.data.frame(do.call("rbind",apply(both[,o.ind],2,vcf.cov.ind)))
#o.icov$Method<-rep("oRAD",nrow(o.icov))
#d.icov<-as.data.frame(do.call("rbind",apply(both[,d.ind],2,vcf.cov.ind)))
#d.icov$Method<-rep("dRAD",nrow(d.icov))
#ind.cov<-as.data.frame(rbind(o.icov,d.icov))
#ind.cov$NumReads<-as.numeric(ind.cov$AvgCovTot)*as.numeric(ind.cov$NumPresent)
#ic<-as.data.frame(apply(ind.cov,2,as.numeric))
#rownames(ic)<-rownames(ind.cov)
#ic$Method<-as.factor(ind.cov$Method)
##Coverage by individual--assembled separately
#orad.icov<-as.data.frame(do.call("rbind",apply(orad[,10:ncol(orad)],2,vcf.cov.ind)))
#oic<-apply(orad.icov,2,as.numeric)
#rownames(oic)<-rownames(orad.icov)
#drad.icov<-as.data.frame(do.call("rbind",apply(drad[,10:ncol(drad)],2,vcf.cov.ind)))
#dic<-apply(drad.icov,2,as.numeric)
#rownames(dic)<-rownames(drad.icov)
##write these to file so I can easily pick back up later
#write.table(oic,"orad_coverage.csv",quote=T,row.names=T,col.names=T,sep='\t')
#write.table(dic,"drad_coverage.csv",quote=T,row.names=T,col.names=T,sep='\t')
#write.table(ic,"both_coverage.csv",quote=T,row.names=T,col.name=T,sep='\t')
oic<-read.delim("orad_coverage.csv")
dic<-read.delim("drad_coverage.csv")
ic<-read.delim("both_coverage.csv")
#subset
o.icov<-as.data.frame(do.call("rbind",apply(both.sub[,o.ind],2,vcf.cov.ind)))
o.icov$Method<-rep("oRAD",nrow(o.icov))
d.icov<-as.data.frame(do.call("rbind",apply(both.sub[,d.ind],2,vcf.cov.ind)))
d.icov$Method<-rep("dRAD",nrow(d.icov))
ind.cov<-as.data.frame(rbind(o.icov,d.icov))
ind.cov$NumReads<-as.numeric(ind.cov$AvgCovTot)*as.numeric(ind.cov$NumPresent)
ic<-as.data.frame(apply(ind.cov,2,as.numeric))
rownames(ic)<-rownames(ind.cov)
ic$Method<-as.factor(ind.cov$Method)
#Coverage by individual--assembled separately
orad.icov<-as.data.frame(do.call("rbind",apply(orad.sub[,10:ncol(orad.sub)],2,vcf.cov.ind)))
oic<-apply(orad.icov,2,as.numeric)
rownames(oic)<-rownames(orad.icov)
drad.icov<-as.data.frame(do.call("rbind",apply(drad.sub[,10:ncol(drad.sub)],2,vcf.cov.ind)))
dic<-apply(drad.icov,2,as.numeric)
rownames(dic)<-rownames(drad.icov)
#write these to file so I can easily pick back up later
write.table(oic,"orad_coverage_subset.csv",quote=T,row.names=T,col.names=T,sep='\t')
write.table(dic,"drad_coverage_subset.csv",quote=T,row.names=T,col.names=T,sep='\t')
write.table(ic,"both_coverage_subset.csv",quote=T,row.names=T,col.name=T,sep='\t')

oic<-read.csv("orad_coverage_subset.csv",row.names=1,header=T,sep='\t')
dic<-read.csv("drad_coverage_subset.csv",row.names=1,header=T,sep='\t')
ic<-read.csv("both_coverage_subset.csv",row.names=1,header=T,sep='\t')

ic.o<-ic[rownames(ic) %in% rownames(oic),]
ic.d<-ic[rownames(ic) %in% rownames(dic),]

mo<-merge(ic.o,oic,by=0)
md<-merge(ic.d,dic,by=0)
oc<-data.frame(NumReadsTogether=mo$NumReads.x,NumReadsAlone=mo$NumReads.y,row.names=mo$Row.names)
dc<-data.frame(NumReadsTogether=md$NumReads.x,NumReadsAlone=md$NumReads.y,row.names=md$Row.names)
t.test(oc$NumReadsTogether,oc$NumReadsAlone,paired=T,alternative="less")
t.test(dc$NumReadsTogether,dc$NumReadsAlone,paired=T,alternative="greater")
wilcox.test(oc$NumReadsTogether,oc$NumReadsAlone,paired=T,alternative="less")
wilcox.test(dc$NumReadsTogether,dc$NumReadsAlone,paired=T,alternative="greater")
oic<-as.data.frame(oic)
dic<-as.data.frame(dic)
oic<-oic[o.ind,]#they had SNP as a row
dic<-dic[d.ind,]
ind.cov<-data.frame(LibraryPrep=c(rep("sdRAD",nrow(oic)),rep("ddRAD",nrow(dic)),
    rep("sdRAD",nrow(oic)),rep("ddRAD",nrow(dic))), 
  Assembly=c(rep("Alone",nrow(oic)),rep("Alone",nrow(dic)),
    rep("Together",nrow(ic))),
  NumReads=c(oic$NumReads,dic$NumReads,ic$NumReads),
  NumMissing=c(oic$NumMissing,dic$NumMissing,ic$NumMissing),
  NumPresent=c(oic$NumPresent,dic$NumPresent,ic$NumPresent))

summary(lm(log(ind.cov$NumReads)~ind.cov$LibraryPrep*ind.cov$Assembly))
summary(glm(log(ind.cov$NumReads)~ind.cov$LibraryPrep*ind.cov$Assembly))
summary(glm(cbind(ind.cov$NumMissing,ind.cov$NumPresent)~ind.cov$LibraryPrep*ind.cov$Assembly,family=binomial))
summary(aov(log(ind.cov$NumReads)~ind.cov$LibraryPrep*ind.cov$Assembly))
TukeyHSD(aov(log(ind.cov$NumReads)~ind.cov$LibraryPrep*ind.cov$Assembly))

####PLOT: Fig 1. Coverage Assembly Method Comp####
jpeg("CoverageAssemblyMethodComp_Subset_rev.jpeg",height=10.5,width=8,units="in",res=300)
par(mfrow=c(3,2),oma=c(1,1,1,1),mar=c(2,2,2,2))
par(mar=c(2,3,2,1))
hist(oc$NumReadsTogether, col=alpha(sdtog.col,0.5), ylim=c(0,150), breaks=seq(0,3500000,250000),
     main="sdRAD-seq",axes=F,xlim=c(0,3500000),border=alpha(sdtog.col,0.5))
hist(oc$NumReadsAlone,col=alpha(sdsep.col,0.5), add=T,breaks=seq(0,3500000,250000),density=20)
legend("topright",c("Assembled Together", "Assembled Separately"),
       fill=c(alpha(sdtog.col,0.5),alpha(sdsep.col,0.5)),bty='n',
       density=c(NA,20),border=c(alpha(sdtog.col,0.5),alpha(sdsep.col,0.5)))
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Number of Individuals",2.5,outer=F,line=2,cex=0.75)
par(mar=c(2,2,2,2))
hist(dc$NumReadsTogether, col=alpha(ddtog.col,0.5), ylim=c(0,150), 
     breaks=seq(0,1500000,100000),border=alpha(ddtog.col,0.5),
     main="ddRAD-seq",axes=F)
hist(dc$NumReadsAlone,col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,1500000,100000),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
legend("topright",c("Assembled Together", "Assembled Separately"),
       fill=c(alpha(ddtog.col,0.5),alpha(ddsep.col,0.5)),bty='n',density = c(NA,20),
       border=c(alpha(ddtog.col,0.5),alpha(ddsep.col,0.5)))
mtext("Total Number of Reads",1,outer=T,cex=0.75,line=-52)#-29 for R

par(mar=c(2,3,2,1))
hist(bo.cov$AvgCovTotal,col=alpha(sdtog.col,0.5),axes=F,xlab="",ylab="",main="",
     xlim=c(0,20080),ylim=c(0,150000),breaks=seq(0,21000,1000),border=alpha(sdtog.col,0.5))
hist(o.cov$AvgCovTotal, col=alpha(sdsep.col,0.5), add=T,breaks=seq(0,21000,1000),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Number of Loci",2,outer=F,line=2.5,cex=0.75)
par(mar=c(2,2,2,2))
hist(bd.cov$AvgCovTotal, col=alpha(ddtog.col,0.5),main="",axes=F,ylab="",xlab="",
     xlim=c(0,2450),ylim=c(0,150000),breaks=seq(0,2600,150),border=alpha(ddtog.col,0.5))
hist(d.cov$AvgCovTotal, col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,2600,150),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Average Number of Reads Per Individual Per SNP",1,outer=T,line=-27,cex=0.75)#-15 for R

par(mar=c(2,3,2,1))
hist(bo.cov$AvgCovRatio,col=alpha(sdtog.col,0.5),main="",axes=F,
     xlab="",ylab="",breaks=seq(0,600,35),ylim=c(0,150000),xlim=c(0,565),border=alpha(sdtog.col,0.5))
hist(o.cov$AvgCovRatio, col=alpha(sdsep.col,0.5), add=T,breaks=seq(0,600,35),density=20)
axis(1,pos=0)
axis(2,pos=-5,las=1)
mtext("Number of Loci",2,outer=F,line=2.5,cex=0.75)
par(mar=c(2,2,2,2))
hist(bd.cov$AvgCovRatio, col=alpha(ddtog.col,0.5),main="",axes=F,
     ylab="",xlab="",breaks=seq(0,2000,125),ylim=c(0,150000),xlim=c(0,2000),border=alpha(ddtog.col,0.5))
hist(d.cov$AvgCovRatio, col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,2000,125),density=20)
axis(1,pos=0)
axis(2,pos=-2,las=1)
mtext("Average Number of Reads in Ref/Avg Number of Reads in Alt",1,outer=T,cex=0.75)
dev.off()

###
jpeg("CoverageAssemblyMethodComp_Subset_rev_small.jpeg",height=10.5,width=8,units="in",res=300)
par(mfrow=c(3,2),oma=c(1,1,1,1),mar=c(2,2,2,2))
par(mar=c(2,3,2,1))
hist(oc$NumReadsTogether, col=alpha(sdtog.col,0.5), ylim=c(0,50),
     breaks=seq(0,3500000,25000),border=alpha(sdtog.col,0.5),
     main="sdRAD-seq",axes=F,xlim=c(0,3500000))
hist(oc$NumReadsAlone,colalpha(sdsep.col,0.5), add=T,breaks=seq(0,3500000,25000),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Number of Individuals",2.5,outer=F,line=2,cex=0.75)
par(mar=c(2,2,2,2))
hist(dc$NumReadsTogether, col=alpha(ddtog.col,0.5), ylim=c(0,50), 
     breaks=seq(0,1500000,10000),border=alpha(ddtog.col,0.5),
     main="ddRAD-seq",axes=F)
hist(dc$NumReadsAlone,col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,1500000,10000),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
legend("topright",c("Assembled Together", "Assembled Separately"),
       fill=c(alpha(sdtog.col,0.5),alpha(sdsep.col,0.5)),bty='n',
       density=c(NA,20),border=c(alpha(sdtog.col,0.5),alpha(sdsep.col,0.5)))
mtext("Total Number of Reads",1,outer=T,cex=0.75,line=-52)#-29 for R

par(mar=c(2,3,2,1))
hist(bo.cov$AvgCovTotal,col=alpha(sdtog.col,0.5),axes=F,xlab="",ylab="",main="",
     xlim=c(0,20080),ylim=c(0,150000),breaks=seq(0,21000,100),border=alpha(sdtog.col,0.5))
hist(o.cov$AvgCovTotal, col=alpha(sdsep.col,0.5), add=T,breaks=seq(0,21000,100),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Number of Loci",2,outer=F,line=2.5,cex=0.75)
par(mar=c(2,2,2,2))
hist(bd.cov$AvgCovTotal, col=alpha(ddtog.col,0.5),main="",axes=F,ylab="",xlab="",
     xlim=c(0,2450),ylim=c(0,150000),breaks=seq(0,2600,15),border=alpha(ddtog.col,0.5))
hist(d.cov$AvgCovTotal, col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,2600,15),density=20)
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Average Number of Reads Per Individual Per SNP",1,outer=T,line=-27,cex=0.75)#-15 for R

par(mar=c(2,3,2,1))
hist(bo.cov$AvgCovRatio,col=alpha(sdtog.col,0.5),main="",axes=F,
     xlab="",ylab="",breaks=seq(0,600,3.5),ylim=c(0,50000),xlim=c(0,565),border=alpha(sdtog.col,0.5))
hist(o.cov$AvgCovRatio, col=alpha(sdsep.col,0.5), add=T,breaks=seq(0,600,3.5),density=20)
axis(1,pos=0)
axis(2,pos=-5,las=1)
mtext("Number of Loci",2,outer=F,line=2.5,cex=0.75)
par(mar=c(2,2,2,2))
hist(bd.cov$AvgCovRatio, col=alpha(ddtog.col,0.5),main="",axes=F,
     ylab="",xlab="",breaks=seq(0,2000,12.5),ylim=c(0,50000),xlim=c(0,2000),border=alpha(ddtog.col,0.5))
hist(d.cov$AvgCovRatio, col=alpha(ddsep.col,0.5), add=T,breaks=seq(0,2000,12.5),density=20)
axis(1,pos=0)
axis(2,pos=-2,las=1)
mtext("Average Number of Reads in Ref/Avg Number of Reads in Alt",1,outer=T,cex=0.75)
dev.off()
####Variance in coverage####
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

##The real analysis
lcv.comp<-data.frame(LibraryPrep=c(rep("sdRAD",nrow(o.cov)),rep("ddRAD",nrow(d.cov)),
                                  rep("sdRAD",nrow(bo.cov)),rep("ddRAD",nrow(bd.cov))), 
                    Assembly=c(rep("Alone",nrow(o.cov)),rep("Alone",nrow(d.cov)),
                               rep("Together",nrow(bo.cov)),rep("Together",nrow(bd.cov))),
                    CovVariance=c(o.cov$CovVariance,d.cov$CovVariance,bo.cov$CovVariance,bd.cov$CovVariance),
                    PropHet=c(o.cov$PropHet,d.cov$PropHet,bo.cov$PropHet,bd.cov$PropHet))

summary(aov(log(lcv.comp$CovVariance+1)~lcv.comp$LibraryPrep*lcv.comp$Assembly))
TukeyHSD(aov(log(lcv.comp$CovVariance+1)~lcv.comp$LibraryPrep*lcv.comp$Assembly))
summary(aov(lcv.comp$PropHet~lcv.comp$LibraryPrep*lcv.comp$Assembly))
TukeyHSD(aov(lcv.comp$PropHet~lcv.comp$LibraryPrep*lcv.comp$Assembly))

####PLOT: Fig 2. Variance in Coverage####
jpeg("VarianceInCov.jpeg",height=8.5,width=7,units="in",res=300)
par(mfrow=c(2,1),oma=c(1,2,1,1),mar=c(2,2,2,2))
boxplot(log(lcv.comp$CovVariance+1)~lcv.comp$LibraryPrep*lcv.comp$Assembly,col=c(ddsep.col,sdsep.col,ddtog.col,sdtog.col),
        names=F,xaxt='n',ylim=c(0,18))
axis(1,at=c(1.5,3.5),c("Alone","Together"))
mtext("ln(Variance in Coverage)",2,outer=F,line=2)
legend("top",ncol=2,c("ddRAD Alone","sdRAD Alone","ddRAD Together", "sdRAD Together"),
       pch=15,col=c(ddsep.col,sdsep.col,ddtog.col,sdtog.col),bty='n')
boxplot(lcv.comp$PropHet~lcv.comp$LibraryPrep*lcv.comp$Assembly,col=c(ddsep.col,sdsep.col,ddtog.col,sdtog.col),
        names=F,xaxt='n')
mtext("Proportion Heterozygotes",2,outer=F,line=2)
axis(1,at=c(1.5,3.5),c("Alone","Together"))
dev.off()

###########################LOOKING FOR THE SAME LOCI#################################
#sdRAD-ddRAD
od.loci<-merge(orad[,locus.info],drad[,locus.info],"SNP")
#od.loci<-merge(orad.sub[,locus.info],drad.sub[,locus.info],"SNP")
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
#ob.loci<-merge(orad.sub[,locus.info],both.sub[,locus.info],"SNP")#subset
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
#db.loci<-merge(drad.sub[,locus.info],both.sub[,locus.info],"SNP")#subset
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
#sdb.loci<-merge(od.loci,both.sub[,locus.info],"SNP")#subset
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
#-using the subsetted files from choose.one.snp
o.share<-orad[orad$SNP %in% od.loci$SNP,]
d.share<-drad[drad$SNP %in% od.loci$SNP,]
od.vcf<-merge(orad,drad,"SNP")
#subset
#o.share<-orad.sub[orad.sub$SNP %in% od.loci$SNP,]
#d.share<-drad.sub[drad.sub$SNP %in% od.loci$SNP,]
#od.sub<-merge(orad.sub,drad.sub,"SNP")
#these are already made
#od.fst<-do.call("rbind",apply(o.share,1,fst.two.vcf,vcf2=d.share,match.index="SNP",cov.thresh=0.5))
#od.fst$SNP<-paste(od.fst$Chrom,as.numeric(as.character(od.fst$Pos)),sep=".")
#od.fst<-do.call("rbind",apply(o.share,1,fst.two.vcf,vcf2=d.share,match.index="SNP"))
#write.table(od.fst,"orad-drad.fst.txt",col.names=T,row.names=F,quote=F,sep='\t')
od.fst<-read.table("orad-drad.fst.txt",header=T,sep='\t')
o.cov.pass<-o.cov[o.cov$AvgCovTotal > 3 & o.cov$AvgCovTotal <= 50, "SNP"]
d.cov.pass<-d.cov[d.cov$AvgCovTotal > 3 & d.cov$AvgCovTotal <= 50, "SNP"]
sep.cov.pass<-o.cov.pass[o.cov.pass %in% d.cov.pass]
summary(od.fst[od.fst$SNP %in% sep.cov.pass,"Fst"])

#FROM SAME ANALYSIS-using the subsetted files from choose.one.snp
both$SNP<-paste(both$`#CHROM`,both$POS,sep=".")
both.o<-cbind(both[,locus.info],both[,o.ind])
both.d<-cbind(both[,locus.info],both[,d.ind])
od.both.fst<-do.call("rbind",apply(both.o,1,fst.two.vcf,vcf2=both.d,match.index="SNP",cov.thresh=0.5))
od.both.fst$SNP<-paste(od.both.fst$Chrom,as.numeric(as.character(od.both.fst$Pos)),sep=".")
od.fst.1<-od.both.fst[od.both.fst$Fst ==1,]
bd.cov$SNP<-paste(bd.cov$Chrom,bd.cov$Pos,sep=".")
bo.cov$SNP<-paste(bo.cov$Chrom,bo.cov$Pos,sep=".")
bd.cov.pass<-bd.cov[bd.cov$AvgCovTotal > 3 & bd.cov$AvgCovTotal <=50,"SNP"]
bo.cov.pass<-bo.cov[bo.cov$AvgCovTotal > 3 & bo.cov$AvgCovTotal <=50,"SNP"]
cov.pass<-bd.cov.pass[bd.cov.pass %in% bo.cov.pass]
summary(od.both.fst[od.both.fst$SNP %in% cov.pass,"Fst"])
#sub.fsts.both<-do.call("rbind",apply(both,1,fst.one.vcf,group1=c(locus.info,o.ind),group2=c(locus.info,d.ind),cov.thresh=0.5))
#sub.fsts.both$SNP<-paste(sub.fsts.both$Chrom,as.numeric(as.character(sub.fsts.both$Pos)),sep=".")
#write.table(sub.fsts.both,"sub.fsts.both.txt",sep='\t',col.names=T,row.names=F,quote=F)

####PLOT: Fig 3. sd-ddRAD Fsts####
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
scaffs<-levels(as.factor(od.fst[,"Chrom"]))
scaffs[1:22]<-lgs

jpeg("sd-dd_Fst_subset.jpeg",height=10,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(2,2,2,2),oma=c(1,2,1,1))
##Separate
od<-fst.plot(od.fst[!is.na(od.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
  fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
  groups=as.factor(scaffs[scaffs %in% levels(factor(od.fst$Chrom[!is.na(od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(od[od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(od[od$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately",3,cex=0.75)
odc<-fst.plot(od.fst[!is.na(od.fst$Fst) & od.fst$SNP %in% cov.pass,],ci.dat=c(-10,10),sig.col=c("black","black"), 
             fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
             groups=as.factor(scaffs[scaffs %in% levels(factor(od.fst$Chrom[!is.na(od.fst$Fst)& od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odc[odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odc[odc$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately, Coverage Filter",3,cex=0.75)

##Together
odb<-fst.plot(od.both.fst[!is.na(od.both.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
              groups=as.factor(scaffs[scaffs %in% levels(factor(od.both.fst[!is.na(od.both.fst$Fst),"Chrom"]))]))
odb$Pos<-as.numeric(as.character(odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odb[odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odb[odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together",3,cex=0.75)
odbc<-fst.plot(od.both.fst[!is.na(od.both.fst$Fst) & od.both.fst$SNP %in% sep.cov.pass,],
  ci.dat=c(-10,10),sig.col=c("black","black"), y.lim=c(0,1),
  fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
  groups=as.factor(scaffs[scaffs %in% 
    levels(factor(od.both.fst[!is.na(od.both.fst$Fst) & od.both.fst$SNP %in% sep.cov.pass,"Chrom"]))]))
odbc$Pos<-as.numeric(as.character(odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together, Coverage Filter",3,cex=0.75)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

gwsca<-read.delim("gwsca_fsts_both.txt")

drad.test<-sample(d.ind,120)
drad.tes.fst<-do.call("rbind", apply(both, 1, fst.one.vcf,group1=drad.test[1:60],group2=drad.test[61:120],cov.thresh=0.2))

##############################60 sdRAD and 60 ddRAD##################################
#d.ind.sub<-sample(d.ind,60,replace=F)
#d.share.sub<-d.share[,colnames(d.share) %in% c(locus.info,d.ind.sub)]
#write.table(d.share,"drad.60ind.vcf",sep='\t',quote=F,row.names=F,col.names=T)
#sub.fsts.both<-do.call("rbind",apply(both.sub,1,fst.one.vcf,group1=c(locus.info,o.ind),group2=colnames(d.share.sub),cov.thresh=0.5))
#sub.fsts.both$SNP<-paste(sub.fsts.both$Chrom,as.numeric(as.character(sub.fsts.both$Pos)),sep=".")
#write.table(sub.fsts.both,"sub.fsts.both.txt",sep='\t',col.names=T,row.names=F,quote=F)
#sub.od.fst<-do.call("rbind",apply(o.share,1,fst.two.vcf,vcf2=d.share,match.index="SNP",cov.thresh=0.5))
#sub.od.fst$SNP<-paste(sub.od.fst$Chrom,as.numeric(as.character(sub.od.fst$Pos)),sep=".")
#write.table(sub.od.fst,"sub.od.fst.txt",sep='\t',col.names=T,row.names=F,quote=F)
d.share.sub<-read.table("d.share.sub",sep='\t',header=T)
sub.od.fst<-read.table("sub.od.fst.txt",sep='\t',header=T)
sub.fsts.both<-read.table("sub.fsts.both.txt",sep='\t',header=T)

png("fsts_60ddrad60sdrad.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(2,2,2,2),oma=c(1,2,1,1))
##Separate
od<-fst.plot(sub.od.fst[!is.na(sub.od.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
             fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
             groups=as.factor(scaffs[scaffs %in% levels(factor(sub.od.fst$Chrom[!is.na(sub.od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(od[od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(od[od$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately",3,cex=0.75)
odc<-fst.plot(sub.od.fst[!is.na(sub.od.fst$Fst) & sub.od.fst$SNP %in% cov.pass,],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
              groups=as.factor(scaffs[scaffs %in% levels(factor(sub.od.fst$Chrom[!is.na(sub.od.fst$Fst)& sub.od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odc[odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odc[odc$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately, Coverage Filter",3,cex=0.75)

##Together
odb<-fst.plot(sub.fsts.both[!is.na(sub.fsts.both$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
              groups=as.factor(scaffs[scaffs %in% levels(factor(sub.fsts.both[!is.na(sub.fsts.both$Fst),"Chrom"]))]))
odb$Pos<-as.numeric(as.character(odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odb[odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odb[odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together",3,cex=0.75)
odbc<-fst.plot(sub.fsts.both[!is.na(sub.fsts.both$Fst) & sub.fsts.both$SNP %in% sep.cov.pass,],
               ci.dat=c(-10,10),sig.col=c("black","black"), y.lim=c(0,1),
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(sub.fsts.both[!is.na(sub.fsts.both$Fst) & sub.fsts.both$SNP %in% sep.cov.pass,"Chrom"]))]))
odbc$Pos<-as.numeric(as.character(odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together, Coverage Filter",3,cex=0.75)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()

#Look at coverage in subset
d60.cov<-do.call("rbind",apply(d.share.sub,1,vcf.cov.loc,subset=d.ind.sub))
wilcox.test(d60.cov$PropHet,o.cov$PropHet)

##############################60 ddRAD and 60 ddRAD##################################
#d.ind.sub1<-sample(d.ind,60,replace=F)
#d.ind.sub2<-sample(d.ind[!(d.ind %in% d.ind.sub1)],60,replace=F)
#dsub.fsts.both<-do.call("rbind",apply(both.sub,1,fst.one.vcf,group1=c(locus.info,d.ind.sub2),group2=c(locus.info,d.ind.sub1),cov.thresh=0.5))
#dsub.fsts.both$SNP<-paste(dsub.fsts.both$Chrom,as.numeric(as.character(dsub.fsts.both$Pos)),sep=".")
#write.table(dsub.fsts.both,"dsub.fsts.both.txt",sep='\t',row.names = F,col.names = T,quote=F)
#d.share.sub1<-d.share[,colnames(d.share) %in% c(locus.info,d.ind.sub1)]
#d.share.sub2<-d.share[,colnames(d.share) %in% c(locus.info,d.ind.sub2)]
#dsub.od.fst<-do.call("rbind",apply(d.share.sub1,1,fst.two.vcf,vcf2=d.share.sub2,match.index="SNP",cov.thresh=0.5))
#dsub.od.fst$SNP<-paste(dsub.od.fst$Chrom,as.numeric(as.character(dsub.od.fst$Pos)),sep=".")
#write.table(d.share.sub1,"drad.60sub1.vcf",row.names=F,col.names=T,quote=F,sep='\t')
#write.table(d.share.sub2,"drad.60sub2.vcf",row.names=F,col.names=T,quote=F,sep='\t')
#write.table(dsub.od.fst,"dsub.od.fst.txt",sep='\t',col.names = T,row.names = F,quote=F)
dsub.od.fst<-read.table("dsub.od.fst.txt",header=T,sep='\t')
dsub.fsts.both<-read.table("dsub.fsts.both.txt",header=T,sep='\t')

png("fsts_60ddrad60sdrad.png",height=10,width=7.5,units="in",res=300)
par(mfrow=c(4,1),mar=c(2,2,2,2),oma=c(1,2,1,1))
##Separate
od<-fst.plot(dsub.od.fst[!is.na(dsub.od.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
             fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
             groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.od.fst$Chrom[!is.na(dsub.od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(od[od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(od[od$Chrom ==lgs[i],"Pos"])
}
mtext("Analyzed Separately",3,cex=0.75)
odc<-fst.plot(dsub.od.fst[!is.na(dsub.od.fst$Fst) & dsub.od.fst$SNP %in% cov.pass,],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
              groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.od.fst$Chrom[!is.na(dsub.od.fst$Fst)& dsub.od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odc[odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odc[odc$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately, Coverage Filter",3,cex=0.75)

##Together
odb<-fst.plot(dsub.fsts.both[!is.na(dsub.fsts.both$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
              groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.fsts.both[!is.na(dsub.fsts.both$Fst),"Chrom"]))]))
odb$Pos<-as.numeric(as.character(odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odb[odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odb[odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together",3,cex=0.75)
odbc<-fst.plot(dsub.fsts.both[!is.na(dsub.fsts.both$Fst) & dsub.fsts.both$SNP %in% sep.cov.pass,],
               ci.dat=c(-10,10),sig.col=c("black","black"), y.lim=c(0,1),
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(dsub.fsts.both[!is.na(dsub.fsts.both$Fst) & dsub.fsts.both$SNP %in% sep.cov.pass,"Chrom"]))]))
odbc$Pos<-as.numeric(as.character(odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odbc[odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
mtext("sdRAD-ddRAD Analyzed Together, Coverage Filter",3,cex=0.75)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
dev.off()


####PLOT: Fig 3NEW. Sample size on Fsts between sd and dd RAD####
png("All3Fsts.png",height=7,width=15,units="in",res=300)
par(mfrow=c(4,3),mar=c(2,1,1,0.01),oma=c(1,3,1,0.5))
##ROW 1: Separate
####All sdRAD vs all ddRAD
a.od<-fst.plot.rect(od.fst[!is.na(od.fst$Fst),],pt.col=ddsdsep.col,
             fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
             groups=as.factor(scaffs[scaffs %in% levels(factor(od.fst$Chrom[!is.na(od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(a.od[a.od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(a.od[a.od$Chrom ==lgs[i],"Pos"])
}   
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(a.od$Fst),4))),bquote(.(nrow(a.od))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("60 sdRAD and 384 ddRAD Individuals",3,cex=0.75,line=0.5)
mtext("Analyzed Separately",2,cex=0.75,line=1.5)

####sdRAD and 60 ddRAD
sd.od<-fst.plot.rect(sub.od.fst[!is.na(sub.od.fst$Fst),], pt.col=ddsdsep.col,
                fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
                groups=as.factor(scaffs[scaffs %in% levels(factor(sub.od.fst$Chrom[!is.na(sub.od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(sd.od[sd.od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(sd.od[sd.od$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(sd.od$Fst),4))),
        bquote(.(nrow(sd.od))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("60 sdRAD and 60 ddRAD Individuals",3,cex=0.75, line = 0.5)

####60 ddRAD and 60 ddRAD
dd.od<-fst.plot.rect(dsub.od.fst[!is.na(dsub.od.fst$Fst),], pt.col = ddsep.col,
                fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),
                groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.od.fst$Chrom[!is.na(dsub.od.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(dd.od[dd.od$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(dd.od[dd.od$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(dd.od$Fst),4))),
        bquote(.(nrow(dd.od))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("60 ddRAD and 60 ddRAD Individuals",3,cex=0.75, line = 0.5)

#####ROW 2
####All sdRAD vs all ddRAD
a.odc<-fst.plot.rect(od.fst[!is.na(od.fst$Fst) & od.fst$SNP %in% cov.pass,], 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddsdsep.col,
              groups=as.factor(scaffs[scaffs %in% levels(factor(od.fst$Chrom[!is.na(od.fst$Fst)& od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(a.odc[a.odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(a.odc[a.odc$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(a.odc$Fst),4))),
        bquote(.(nrow(a.odc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("Analyzed Separately,\nCoverage Filter",2,cex=0.75,line=1.5)
####sdRAD vs 60 ddRAD
sd.odc<-fst.plot.rect(sub.od.fst[!is.na(sub.od.fst$Fst) & sub.od.fst$SNP %in% cov.pass,], 
                 fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddsdsep.col,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(sub.od.fst$Chrom[!is.na(sub.od.fst$Fst)& sub.od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(sd.odc[sd.odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(sd.odc[sd.odc$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(sd.odc$Fst),4))),
        bquote(.(nrow(sd.odc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)


####60 ddRAD and 60 ddRAD
dd.odc<-fst.plot.rect(dsub.od.fst[!is.na(dsub.od.fst$Fst) & dsub.od.fst$SNP %in% cov.pass,],
                 fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddsep.col,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.od.fst$Chrom[!is.na(dsub.od.fst$Fst)& dsub.od.fst$SNP %in% cov.pass]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(dd.odc[dd.odc$Chrom ==lgs[i],"Pos"]),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(dd.odc[dd.odc$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(dd.odc$Fst),4))),
        bquote(.(nrow(dd.odc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)


#####ROW 3: Together
####All sdRAD vs all ddRAD
a.odb<-fst.plot.rect(od.both.fst[!is.na(od.both.fst$Fst),],
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddsdtog.col,
              groups=as.factor(scaffs[scaffs %in% levels(factor(od.both.fst[!is.na(od.both.fst$Fst),"Chrom"]))]))
a.odb$Pos<-as.numeric(as.character(a.odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(a.odb[a.odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(a.odb[a.odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(a.odb$Fst),4))),
        bquote(.(nrow(a.odb))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("Analyzed Together",2,line=1.5,cex=0.75)

####All sdRAD vs 60 ddRAD
sd.odb<-fst.plot.rect(sub.fsts.both[!is.na(sub.fsts.both$Fst),],
                 fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddsdtog.col,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(sub.fsts.both[!is.na(sub.fsts.both$Fst),"Chrom"]))]))
sd.odb$Pos<-as.numeric(as.character(sd.odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(sd.odb[sd.odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(sd.odb[sd.odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(sd.odb$Fst),4))),
         bquote(.(nrow(sd.odb))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)

####60 ddRAD and 60 ddRAD
dd.odb<-fst.plot.rect(dsub.fsts.both[!is.na(dsub.fsts.both$Fst),],
                 fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(0,1),pt.col = ddtog.col,
                 groups=as.factor(scaffs[scaffs %in% levels(factor(dsub.fsts.both[!is.na(dsub.fsts.both$Fst),"Chrom"]))]))
dd.odb$Pos<-as.numeric(as.character(dd.odb$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(dd.odb[dd.odb$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(dd.odb[dd.odb$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(dd.odb$Fst),4))),
        bquote(.(nrow(dd.odb))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)

#####ROW 4: Together
####All sdRAD vs all ddRAD
a.odbc<-fst.plot.rect(od.both.fst[!is.na(od.both.fst$Fst) & od.both.fst$SNP %in% sep.cov.pass,],
              y.lim=c(0,1),pt.col = ddsdtog.col,
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(od.both.fst[!is.na(od.both.fst$Fst) & od.both.fst$SNP %in% sep.cov.pass,"Chrom"]))]))
a.odbc$Pos<-as.numeric(as.character(a.odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(a.odbc[a.odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(a.odbc[a.odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(a.odbc$Fst),4))),
        bquote(.(nrow(a.odbc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
mtext("Analyzed Together,\nCoverage Filter",2,line=1.5,cex=0.75)

####All sdRAD vs 60 ddRAD
sd.odbc<-fst.plot.rect(sub.fsts.both[!is.na(sub.fsts.both$Fst) & sub.fsts.both$SNP %in% sep.cov.pass,],
               y.lim=c(0,1),pt.col = ddsdtog.col,
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(sub.fsts.both[!is.na(sub.fsts.both$Fst) & sub.fsts.both$SNP %in% sep.cov.pass,"Chrom"]))]))
sd.odbc$Pos<-as.numeric(as.character(sd.odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(sd.odbc[sd.odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(sd.odbc[sd.odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(sd.odbc$Fst),4))),
        bquote(.(nrow(sd.odbc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)

###60 ddRAD vs 60 ddRAD
dd.odbc<-fst.plot.rect(dsub.fsts.both[!is.na(dsub.fsts.both$Fst) & dsub.fsts.both$SNP %in% sep.cov.pass,],
              y.lim=c(0,1),pt.col = ddtog.col,
               fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,
               groups=as.factor(scaffs[scaffs %in% 
                                         levels(factor(dsub.fsts.both[!is.na(dsub.fsts.both$Fst) & dsub.fsts.both$SNP %in% sep.cov.pass,"Chrom"]))]))
dd.odbc$Pos<-as.numeric(as.character(dd.odbc$Pos))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(dd.odbc[dd.odbc$Chrom ==lgs[i],"Pos"],na.rm=T),y=-0.05,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(dd.odbc[dd.odbc$Chrom ==lgs[i],"Pos"],na.rm = T)
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(dd.odbc$Fst),4))),
        bquote(.(nrow(dd.odbc))~" SNPs"))
legend("top",legend=as.expression(lgnd),bg=alpha("grey94",0.3),box.lty=0)
legend(x=2500,y=0.75,bg=alpha("grey94",0.5),pch=19,col=c(ddsdsep.col,ddsdtog.col,ddsep.col,ddtog.col),box.lty=0,
       c("sdRAD-ddRAD Analyzed Separately","sdRAD-ddRAD Analyzed Together","ddSNPs Analyzed Separately","ddSNPs Analyzed Together"),
       ncol=1)
mtext(expression(italic(F)[ST]),2,outer=T,cex=0.75)
mtext("Linkage Group",1,outer=T,cex=0.75)
dev.off()


fst.comp<-data.frame(Filtering=c(rep("Unfiltered",nrow(a.od)),rep("Filtered",nrow(a.odc)),
                                 rep("Unfiltered",nrow(a.odb)),rep("Filtered",nrow(a.odbc))), 
                     Assembly=c(rep("Alone",nrow(a.od)),rep("Alone",nrow(a.odc)),
                                rep("Together",nrow(a.odb)),rep("Together",nrow(a.odbc))),
                     Fsts=c(a.od$Fst,a.odc$Fst,a.odb$Fst,a.odbc$Fst))
fst.comp.aov<-aov(Fsts~Filtering*Assembly,dat=fst.comp)
summary(fst.comp.aov)
TukeyHSD(fst.comp.aov)


fst.aov.dat<-data.frame(Fst=c(a.od$Fst,sd.od$Fst,dd.od$Fst,a.odb$Fst,sd.odb$Fst,dd.odb$Fst,a.odc$Fst,
                              sd.odc$Fst,dd.odc$Fst,a.odbc$Fst,sd.odbc$Fst,dd.odbc$Fst),
                        SampleSize=c(rep("All",nrow(a.od)),rep("Each60",nrow(sd.od)),rep("Each60",nrow(dd.od)),
                               rep("All",nrow(a.odb)),rep("Each60",nrow(sd.odb)),rep("Each60",nrow(dd.odb)),
                               rep("All",nrow(a.odc)),rep("Each60",nrow(sd.odc)),rep("Each60",nrow(dd.odc)),
                               rep("All",nrow(a.odbc)),rep("Each60",nrow(sd.odbc)),rep("Each60",nrow(dd.odbc))),
                        Comparison=c(rep("sd-dd",nrow(a.od)),rep("sd-dd",nrow(sd.od)),rep("dd-dd",nrow(dd.od)),
                                     rep("sd-dd",nrow(a.odb)),rep("sd-dd",nrow(sd.odb)),rep("dd-dd",nrow(dd.odb)),
                                     rep("sd-dd",nrow(a.odc)),rep("sd-dd",nrow(sd.odc)),rep("dd-dd",nrow(dd.odc)),
                                     rep("sd-dd",nrow(a.odbc)),rep("sd-dd",nrow(sd.odbc)),rep("dd-dd",nrow(dd.odbc))),
                        Analysis=c(rep("UnfilteredSeparate",(nrow(a.od)+nrow(sd.od)+nrow(dd.od))),
                                   rep("UnfilteredTogether",(nrow(a.odb)+nrow(sd.odb)+nrow(dd.odb))),
                                   rep("FilteredSeparate",(nrow(a.odc)+nrow(sd.odc)+nrow(dd.odc))),
                                   rep("FilteredTogether",(nrow(a.odbc)+nrow(sd.odbc)+nrow(dd.odbc)))))
fst.ss.aov<-aov(Fst~SampleSize*Analysis,dat=fst.aov.dat[fst.aov.dat$Comparison == "sd-dd",])
fst.cp.aov<-aov(Fst~Comparison*Analysis,dat=fst.aov.dat)

##############################SCA##################################
#BOTH
setwd("both_sca/")
both.sexsel<-read.table("both.sexsel.txt",header=T,sep='\t')
both.viasel<-read.table("both.viasel.txt",header=T,sep='\t')
drad.sexsel<-read.table("drad.sexsel.txt",header=T,sep='\t')
drad.viasel<-read.table("drae.viasel.txt",header=T,sep='\t')
orad.viasel<-read.table("orad.viasel.txt",header=T,sep='\t')
#INFER MATERNAL ALLELES

#SCA
locus.info<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
vcf1<-parse.vcf("both_maternal.vcf")
#rename moms
names1<-colnames(vcf1)
names1[grep("MOM",names1)]<-gsub("MOM.*PRM(\\d{3})(_align)?","MOM\\1",names1[grep("MOM",names1)])
colnames(vcf1)<-names1
vcf2.gt<-extract.gt.vcf(both)
merge<-merge.vcfs(vcf1,both)
both.keep<-b.cov[b.cov$AvgCovTotal > 5 & b.cov$AvgCovTotal <= 20,"SNP"]
merge$SNP<-paste(as.character(merge$CHROM),as.character(merge$POS),sep=".")
merge<-merge[merge$SNP %in% both.keep,]
females<-colnames(merge[grep("FEM",colnames(merge))])
males<-c(colnames(merge[grep("PRM",colnames(merge))]),colnames(merge[grep("NPM",colnames(merge))]))
moms<-colnames(merge[grep("MOM",colnames(merge))])
#both.sexsel<-gwsca(merge,locus.info,females,moms)
#both.sexsel$index<-paste(both.sexsel$Chrom,both.sexsel$Pos,sep=".")
bss.sig<-both.sexsel$index[both.sexsel$Chi.p.adj <= 0.05]
#both.viasel<-gwsca(merge,locus.info,females,males)
#both.viasel$index<-paste(both.viasel$Chrom,both.viasel$Pos,sep=".")
bvs.sig<-both.viasel$index[both.viasel$Chi.p.adj <= 0.05]
#write.table(both.sexsel,"both.sexsel.txt",col.names=T,row.names=F,quote=F,sep='\t')
#write.table(both.viasel,"both.viasel.txt",col.names=T,row.names=F,quote=F,sep='\t')

#DRAD
drad.merge<-read.table("../biallelic/biallelic_merge.vcf",sep='\t',header=T)
drad.keep<-d.cov[d.cov$AvgCovTotal > 5 & d.cov$AvgCovTotal <= 20,"SNP"]
drad.merge$SNP<-paste(drad.merge$CHROM,drad.merge$POS,sep=".")
drad.merge<-drad.merge[drad.merge$SNP %in% drad.keep,]
#rename moms
#names1<-colnames(drad.merge)
#names1[grep("MOM",names1)]<-gsub("MOM.*PRM(\\d{3})(_align)?","MOM\\1",names1[grep("MOM",names1)])
#colnames(drad.merge)<-names1
#locus.info<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
#females<-colnames(drad.merge[grep("FEM",colnames(drad.merge))])
#males<-c(colnames(drad.merge[grep("PRM",colnames(drad.merge))]),colnames(drad.merge[grep("NPM",colnames(drad.merge))]))
#moms<-colnames(drad.merge[grep("MOM",colnames(drad.merge))])
#drad.sexsel<-gwsca(drad.merge,locus.info,females, moms)
#drad.sexsel$index<-paste(drad.sexsel$Chrom,drad.sexsel$Pos,sep=".")
dss.sig<-drad.sexsel$index[drad.sexsel$Chi.p.adj <= 0.05]
#drad.viasel<-gwsca(drad.merge,locus.info,females,males)
#drad.viasel$index<-paste(drad.viasel$Chrom,drad.viasel$Pos,sep=".")
dvs.sig<-drad.viasel$index[drad.viasel$Chi.p.adj <= 0.05]
#write.table(drad.sexsel,"drad.sexsel.txt",col.names=T,row.names=F,quote=F,sep='\t')
#write.table(drad.viasel,"drae.viasel.txt",col.names=T,row.names=F,quote=F,sep='\t')

#SRAD
locus.info<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
orad.keep<-o.cov[o.cov$AvgCovTotal > 5 & o.cov$AvgCovTotal <= 20,"SNP"]
orad<-orad[orad$SNP %in% orad.keep,]
females<-colnames(orad[grep("FEM",colnames(orad))])
males<-colnames(orad[grep("PRM",colnames(orad))])
orad.viasel<-gwsca(orad,locus.info,females,males)
#orad.viasel$index<-paste(orad.viasel$Chrom,orad.viasel$Pos,sep=".")
ovs.sig<-orad.viasel$index[orad.viasel$Chi.p.adj <= 0.05]
#write.table(orad.viasel,"orad.viasel.txt",col.names=T,row.names=F,quote=F,sep='\t')

#STATISTICAL ANALYSIS
viasel.dat<-data.frame(Fst=c(orad.viasel$Fst,drad.viasel$Fst,both.viasel$Fst),
                       Method=c(rep("sdRAD",nrow(orad.viasel)),rep("ddRAD",nrow(drad.viasel)),rep("Both",nrow(both.viasel))))
via.aov<-aov(Fst~Method,dat=viasel.dat)
sexsel.dat<-data.frame(Fst=c(drad.sexsel$Fst,both.sexsel$Fst),
                       Method=c(rep("ddRAD",nrow(drad.sexsel)),rep("Both",nrow(both.sexsel))))
sex.aov<-aov(Fst~Method,dat=sexsel.dat)
#sca.dat<-data.frame(Fst=c(orad.viasel$Fst,drad.viasel$Fst,both.viasel$Fst,drad.sexsel$Fst,both.sexsel$Fst),
#  Method=c(rep("sdRAD",nrow(orad.viasel)),rep("ddRAD",nrow(drad.viasel)),rep("Both",nrow(both.viasel)),
#           rep("ddRAD",nrow(drad.sexsel)),rep("Both",nrow(both.sexsel))),
#  Analysis=c(rep("Viability",nrow(orad.viasel)),rep("Viability",nrow(drad.viasel)),rep("Viability",nrow(both.viasel)),
#             rep("Sexual",nrow(drad.sexsel)),rep("Sexual",nrow(both.sexsel))))
#sca.aov<-aov(Fst~Method*Analysis,dat=sca.dat)

#check to make sure differences aren't just due to sample size
locus.info<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")
females<-colnames(drad.merge[grep("FEM",colnames(drad.merge))])
males<-c(colnames(drad.merge[grep("PRM",colnames(drad.merge))]),colnames(drad.merge[grep("NPM",colnames(drad.merge))]))
d30mal<-sample(males,30,replace = F)
d30fem<-sample(females,30,replace = F)
drad.viasel.30<-gwsca(drad.merge,locus.info,d30fem,d30mal)
drad.viasel.30$index<-paste(as.character(drad.viasel.30$Chrom),as.numeric(drad.viasel.30$Pos),sep=".")
dvs30.sig<-drad.viasel.30$index[drad.viasel.30$Chi.p.adj <= 0.05]
fst.plot.rect(fst.dat=drad.viasel.30[!is.na(drad.viasel.30$Fst),],pt.col = ddsep.col, 
         fst.name="Fst", chrom.name="Chrom", bp.name="Pos",y.lim=c(0,0.5),axis.size=1,
         groups=as.factor(scaffs[scaffs %in% levels(factor(drad.viasel.30$Chrom[!is.na(drad.viasel.30$Fst)]))]))
wilcox.test(drad.viasel.30$Fst,orad.viasel$Fst,alternative="less")
wilcox.test(drad.viasel.30$Fst,drad.viasel$Fst,alternative="greater")
####PLOT: Fig 4. SCA####
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
scaffs<-levels(as.factor(both[,"#CHROM"]))
scaffs[1:22]<-lgs

png("SCA_radseq.png",height=7.5,width=10,units="in",res=300)
par(mfrow=c(2,3),mar=c(2,1.5,2,1),oma=c(1,2,1,0.5))
##ROW 1: Males vs Females
####both sdRAD and ddRAD (RAD)
vs<-fst.plot.rect(fst.dat=both.viasel[!is.na(both.viasel$Fst),],pt.col=ddsdtog.col, 
  fst.name="Fst", chrom.name="Chrom", bp.name="Pos",y.lim=c(0,0.5),axis.size=1,
  groups=as.factor(scaffs[scaffs %in% levels(factor(both.viasel$Chrom[!is.na(both.viasel$Fst)]))]))
points(vs$Pos[vs$index %in% bvs.sig], 
       vs$Fst[vs$index %in% bvs.sig],
       col=viasel.col,pch=19,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(vs[vs$Chrom ==lgs[i],"Pos"]),y=-0.01,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(vs[vs$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(vs$Fst),4))),
        bquote(.(nrow(vs))~" SNPs"), bquote(.(length(bvs.sig))~" Significant SNPs"))
legend("top",inset=c(0,0.05),legend=as.expression(lgnd), bg=alpha("grey94",0.5),box.lty=0)
mtext("Both ddRAD-seq and sdRAD-seq",3,cex=0.75,line=0.5)
mtext(expression(Males-Females~italic(F)[ST]),2,line=2,cex=0.75)

####ddRAD-seq only
vsd<-fst.plot.rect(fst.dat=drad.viasel[!is.na(drad.viasel$Fst),],pt.col = ddsep.col, 
              fst.name="Fst", chrom.name="Chrom", bp.name="Pos",y.lim=c(0,0.5),axis.size=1,
              groups=as.factor(scaffs[scaffs %in% levels(factor(drad.viasel$Chrom[!is.na(drad.viasel$Fst)]))]))
points(vsd$Pos[vsd$Chi.p.adj <= 0.05], 
       vsd$Fst[vsd$Chi.p.adj <= 0.05],
       col=viasel.col,pch=19,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(as.numeric(vsd[vsd$Chrom ==lgs[i],"Pos"])),y=-0.01,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(as.numeric(vsd[vsd$Chrom ==lgs[i],"Pos"]))
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(vsd$Fst),4))),
         bquote(.(nrow(vsd))~" SNPs"), bquote(.(length(vsd$Fst[vsd$Chi.p.adj <= 0.05]))~" Significant SNPs"))
legend("top",inset=c(0,0.05),legend=as.expression(lgnd), bg=alpha("grey94",0.5),box.lty=0)
mtext("ddRAD-seq",3,cex=0.75, line = 0.5)

####sdRAD-seq only
vso<-fst.plot.rect(fst.dat=orad.viasel[!is.na(orad.viasel$Fst),],pt.col=sdsep.col, 
               fst.name="Fst", chrom.name="Chrom", bp.name="Pos", axis.size=1,y.lim=c(0,0.5),
                groups=as.factor(scaffs[scaffs %in% levels(factor(orad.viasel$Chrom[!is.na(orad.viasel$Fst)]))]))
points(vso$Pos[vso$index %in% ovs.sig], 
       vso$Fst[vso$index %in% ovs.sig],
       col=viasel.col,pch=19,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(as.numeric(vso[vso$Chrom ==lgs[i],"Pos"])),y=-0.01,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(as.numeric(vso[vso$Chrom ==lgs[i],"Pos"]))
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(vso$Fst),4))),
        bquote(.(nrow(vso))~" SNPs"), bquote(.(length(ovs.sig))~" Significant SNPs"))
legend("top",inset=c(0,0.05),legend=as.expression(lgnd), bg=alpha("grey94",0.5),box.lty=0)
mtext("sdRAD-seq",3,cex=0.75, line = 0.5)

#####ROW 2: sexual selection
####Both sdRAD-seq and ddRAD-seq
ss<-fst.plot.rect(fst.dat=both.sexsel[!is.na(both.sexsel$Fst),],pt.col=ddsdtog.col, 
             fst.name="Fst", chrom.name="Chrom", bp.name="Pos",axis.size=1,y.lim=c(0,0.5),
                groups=as.factor(scaffs[scaffs %in% levels(factor(both.sexsel$Chrom[!is.na(both.sexsel$Fst)]))]))
points(ss$Pos[ss$index %in% bss.sig], 
       ss$Fst[ss$index %in% bss.sig],
       col=sexsel.col,pch=19,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(ss[ss$Chrom ==lgs[i],"Pos"]),y=-0.01,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(ss[ss$Chrom ==lgs[i],"Pos"])
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(ss$Fst),4))),
        bquote(.(nrow(ss))~" SNPs"), bquote(.(length(bss.sig))~" Significant SNPs"))
legend("top",inset=c(0,0.05),legend=as.expression(lgnd), bg=alpha("grey94",0.5),box.lty=0)
mtext(expression(Inferred~Mothers-Females~italic(F)[ST]),2,line=2,cex=0.75)

####ddRAD-seq only
ssd<-fst.plot.rect(fst.dat=drad.sexsel[!is.na(drad.sexsel$Fst),],pt.col=ddsep.col, 
              fst.name="Fst", chrom.name="Chrom", bp.name="Pos",axis.size=1,y.lim=c(0,0.5),
              groups=as.factor(scaffs[scaffs %in% levels(factor(drad.sexsel$Chrom[!is.na(drad.sexsel$Fst)]))]))
points(ssd$Pos[ssd$index %in% dss.sig], 
       ssd$Fst[ssd$index %in% dss.sig],
       col=sexsel.col,pch=19,cex=0.75)
last<-0
for(i in 1:length(lgs)){
  text(x=mean(as.numeric(ssd[ssd$Chrom %in%lgs[i],"Pos"])),y=-0.01,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(as.numeric(ssd[ssd$Chrom ==lgs[i],"Pos"]))
}
lgnd<-c(bquote("Mean "~italic(F)[ST]~"="~.(round(mean(ssd$Fst),4))),
        bquote(.(nrow(ssd))~" SNPs"), bquote(.(length(dss.sig))~" Significant SNPs"))
legend("top",inset=c(0,0.05),legend=as.expression(lgnd), bg=alpha("grey94",0.5),box.lty=0)

####Empty plot for legend
plot(1, type="n", xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1),axes=F,xaxt='n',yaxt='n')
legend("top",c("sex-biased viability selection","sexual selection",
               "sdRAD and ddRAD Analyzed Together","ddRAD Analyzed Separately","sdRAD Analyzed Separately"),
       col=c(viasel.col,sexsel.col,ddsdtog.col,ddsep.col,sdsep.col),pch=19,
       ncol=1,bty='n')
mtext("Linkage Group",1,outer=T,cex=0.75)
dev.off()

##############################DRAD DIFFERENT PLATES##################################

#d.cov<-do.call("rbind",apply(drad,1,vcf.cov.loc,subset=d.ind))
#d.cov$SNP<-paste(d.cov$Chrom,d.cov$Locus,d.cov$Pos,sep=".")
#d.cov$AvgCovRatio[d.cov$AvgCovRatio=="Inf"]<-0
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

###############################ALLELE FREQUENCIES#####################################
bd.afs<-do.call("rbind",apply(both[,c(locus.info,d.ind)],1,calc.afs.vcf))
bs.afs<-do.call("rbind",apply(both[,c(locus.info,o.ind)],1,calc.afs.vcf))
bd.afs$index<-paste(bd.afs$Chrom,bd.afs$Pos,sep=".")
bs.afs$index<-paste(bs.afs$Chrom,bs.afs$Pos,sep=".")
bod.afs<-merge(bd.afs,bs.afs,by="index")
wilcox.test(bod.afs$RefFreq.x,bod.afs$RefFreq.y,paired=T,alternative="greater")
mean(bs.afs$RefFreq)
mean(bd.afs$RefFreq)

drad.afs<-do.call("rbind",apply(drad,1,calc.afs.vcf))
orad.afs<-do.call("rbind",apply(orad,1,calc.afs.vcf))
drad.afs$index<-paste(drad.afs$Chrom,drad.afs$Pos,sep=".")
orad.afs$index<-paste(orad.afs$Chrom,orad.afs$Pos,sep=".")
od.afs<-merge(drad.afs,orad.afs,by="index")
wilcox.test(od.afs$RefFreq.x,od.afs$RefFreq.y,paired=T,alternative = "less")
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

####PLINK SUBSETS--DON'T USE####
dsub<-read.delim("drad.subset.raw",sep=" ")
dmap<-read.delim("drad.map",header=F)
osub<-read.delim("orad.subset.raw",sep=" ")
omap<-read.delim("orad.map",header=F)
bsub<-read.delim("both.subset.raw",sep=" ")
bmap<-read.delim("both.map",header=F)

#Match dsub and osub
merge.map<-merge(omap,dmap,by="SNP")
merge.map$oLocus<-apply(merge.map,1,function(x){
  x["oLocus"]<-paste("X",x["oLocus"],sep="")
})
merge.map$dLocus<-apply(merge.map,1,function(x){
  x["dLocus"]<-paste("X",x["dLocus"],sep="")
})
#match orad
o.newcol<-colnames(osub)
o.newcol<-unlist(lapply(o.newcol,function(x){
  o.newcol<-gsub("(X\\d+_\\d+)_\\w","\\1",x)
}))
colnames(osub)<-o.newcol
omerge.keep<-merge.map[merge.map$oLocus %in% colnames(osub),]
omerge.keep<-omerge.keep[!duplicated(omerge.keep$oLocus),]
odsub<-osub[,colnames(osub)%in%omerge.keep$oLocus]
odsub<-odsub[,order(colnames(odsub))]
omerge.keep<-omerge.keep[order(omerge.keep$oLocus),]
colnames(odsub)<-omerge.keep[omerge.keep$oLocus %in% colnames(odsub),"SNP"]
mo<-odsub[,!duplicated(colnames(odsub))]
#match drad
d.newcol<-colnames(dsub)
d.newcol<-unlist(lapply(d.newcol,function(x){
  d.newcol<-gsub("(X\\d+_\\d+)_\\w","\\1",x)
}))
colnames(dsub)<-d.newcol
dmerge.keep<-merge.map[merge.map$dLocus %in% colnames(dsub),]
dmerge.keep<-dmerge.keep[!duplicated(dmerge.keep$dLocus),]
ddsub<-dsub[,colnames(dsub)%in%dmerge.keep$dLocus]
ddsub<-ddsub[,order(colnames(ddsub))]
dmerge.keep<-dmerge.keep[order(dmerge.keep$dLocus),]
colnames(ddsub)<-dmerge.keep[dmerge.keep$dLocus %in% colnames(ddsub),"SNP"]
md<-ddsub[,!duplicated(colnames(ddsub))]
mo<-mo[,order(colnames(mo))]
md<-md[,order(colnames(md))]
mergedsub<-data.frame(rbind(cbind(osub[,1:6],mo[,colnames(mo) %in% colnames(md)]),
                            cbind(dsub[,1:6],md[,colnames(md)%in% colnames(mo)])))
#fsts
plink.diff.fst<-fst.one.plink(mergedsub,group1 = osub$IID,group2 = dsub$IID)
plink.both.fst<-fst.one.plink(bsub,group1 = osub$IID,group2 = dsub$IID)
plink.diff.fst$Chrom<-gsub("(\\w+\\d+).\\d+","\\1",plink.diff.fst$Locus)
plink.diff.fst$Pos<-gsub("(\\w+\\d+).(\\d+)","\\2",plink.diff.fst$Locus)

plink.both.fst$Chrom<-apply(plink.both.fst,1,function(x){
  chrom<-bmap[bmap$V2 %in% gsub("X(\\d+_\\d+)_\\w","\\1",x),"V1"]
})
plink.both.fst$Chrom<-factor(plink.both.fst$Chrom)
plink.both.fst$Pos<-apply(plink.both.fst,1,function(x){
  pos<-bmap[bmap$V2 %in% gsub("X(\\d+_\\d+)_\\w","\\1",x),"V4"]
})

#PLOTTING##
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
       "LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
       "LG22")
lgn<-seq(1,22)
scaffs<-levels(as.factor(plink.both.fst[,"Chrom"]))
scaffs[1:22]<-lgs

jpeg("sd-dd_Fst_subset.jpeg",height=10,width=7.5,units="in",res=300)
par(mfrow=c(2,1),mar=c(2,2,2,2),oma=c(1,2,1,1))
od<-fst.plot(plink.both.fst[!is.na(plink.both.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
             fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(-2,0.6),
             groups=as.factor(scaffs[scaffs %in% levels(factor(plink.both.fst$Chrom[!is.na(plink.both.fst$Fst)]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(od[od$Chrom ==lgs[i],"Pos"]),y=-2.1,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(od[od$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Together",3,cex=0.75)
odb<-fst.plot(plink.diff.fst[!is.na(plink.diff.fst$Fst),],ci.dat=c(-10,10),sig.col=c("black","black"), 
              fst.name="Fst",chrom.name="Chrom",bp.name="Pos",axis.size=1,y.lim=c(-1,0.3),
              groups=as.factor(scaffs[scaffs %in% levels(factor(plink.diff.fst[!is.na(plink.diff.fst),"Chrom"]))]))
last<-0
for(i in 1:length(lgs)){
  text(x=mean(odb[odb$Chrom ==lgs[i],"Pos"]),y=-1.1,
       labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=1)
  last<-max(odb[odb$Chrom ==lgs[i],"Pos"])
}
mtext("sdRAD-ddRAD Analyzed Separately",3,cex=0.75)
dev.off()

##STATS
plink.fst<-data.frame(Fst=as.numeric(c(plink.both.fst$Fst,plink.diff.fst$Fst)),
                      Analysis=c(rep("Together",length(plink.both.fst$Fst)),rep("Separately",length(plink.diff.fst$Fst))))

>>>>>>> parent of 70c4bc9... merging
