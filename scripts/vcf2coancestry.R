#Author: Sarah P. Flanagan
#Date: 23 April 2018
#Purpose: to convert a vcf file to Coancestry allele frequency input file format
#I've added the function to gwscaR

setwd("B:/ubuntushare/SCA/")
setwd("results/")

source("../../gwscaR/R/gwscaR.R")
source("../../gwscaR/R/gwscaR_plot.R")
source("../../gwscaR/R/gwscaR_utility.R")
source("../../gwscaR/R/gwscaR_fsts.R")
source("../../gwscaR/R/gwscaR_popgen.R")

vcf<-parse.vcf("drad.sub.vcf")

#calculate missingness at each locus
gts<-extract.gt.vcf(vcf)
prop.miss<-apply(gts[,10:ncol(gts)],1,function(gt){
  n<-length(gt[gt=="./."])
  prop.missing<-n/length(gt)
  return(prop.missing)
})

vcf<-vcf[which(prop.miss<=0.3),]
vcf$`#CHROM`<-gsub("[A-z]+_?(\\d+)","\\1",vcf$`#CHROM`)
colnames(vcf)[10:ncol(vcf)]<-gsub("sample_(\\w\\w)\\w(\\d.*)_align","\\1_\\2",colnames(vcf[10:ncol(vcf)]))

write.table(vcf,"relatedness/drad_miss3.vcf",col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

#Error file
static.error<-data.frame(rep(0.3,nrow(vcf)),rep(0.04,nrow(vcf)),rep(0.2,nrow(vcf)))
write.table(static.error,"relatedness/coancestry_err.txt",sep='\t',col.names=FALSE,row.names = FALSE,quote=FALSE)

#Allele frequency data
vcf2coanAF<-function(vcf,out.name="coancestry_afs.txt"){
  co.afs<-do.call(rbind,apply(vcf,1,function(vcf.row,out.name){
    af<-calc.afs.vcf(vcf.row)
    bases<-c("A","C","G","T")
    rname<-which(bases==vcf.row[["REF"]])
    aname<-which(bases==vcf.row[["ALT"]])
    co.af<-cbind(round(af$RefFreq,4),round(af$AltFreq,4))
    colnames(co.af)<-c(rname,aname)
    suppressWarnings(write.table(co.af,out.name,sep='\t',append = TRUE,quote=FALSE,row.names = FALSE,col.names = TRUE))
    row.names(co.af)<-vcf.row["ID"]
    colnames(co.af)<-NULL
    return(data.frame(co.af))
  },out.name=out.name))
  return(co.afs)
}

co.afs<-vcf2coanAF(vcf,"relatedness/coancestry_afs.txt")

#Genotype data
vcf2coanGT<-function(vcf,out.name="coancestry_gty.txt"){
  gts<-extract.gt.vcf(vcf)
  co.gt<-do.call(cbind,apply(gts,1,function(gt){
    bases<-c("A","C","G","T")
    rname<-which(bases==gt[["REF"]])
    aname<-which(bases==gt[["ALT"]])
    g<-do.call(rbind,lapply(gt[10:length(gt)],function(x) {
      strsplit(as.character(x),"/")[[1]]}))
    #print(rname)
    g[g=="0"]<-rname
    g[g=="1"]<-aname
    g[g=="."]<-0
    return(as.data.frame(g))
  }))
  
  write.table(co.gt,out.name,row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
  return(co.gt)
}

co.gt<-vcf2coanGT(vcf,"relatedness/coancestry_gty.txt")

