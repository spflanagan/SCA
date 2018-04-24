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

out.name<-"coancestry_afs.txt"
#Allele frequency data
vcf2coanAF<-function(vcf,out.name="coancestry_afs.txt"){
  co.afs<-do.call(rbind,apply(vcf,1,function(vcf.row,out.name){
    af<-calc.afs.vcf(vcf.row)
    snp.name<-paste(af$Chrom,trimws(af$Pos),sep=".")
    co.af<-cbind(af$RefFreq,af$AltFreq)
    colnames(co.af)<-c(paste(snp.name,af$Ref,sep="_"),paste(snp.name,af$Alt,sep="_"))
    suppressWarnings(write.table(co.af,out.name,sep='\t',append = TRUE,quote=FALSE,row.names = FALSE,col.names = TRUE))
    row.names(co.af)<-snp.name
    colnames(co.af)<-NULL
    return(data.frame(co.af))
  },out.name=out.name))
  return(co.afs)
}

co.afs<-vcf2coanAF(vcf)
co.afs<-vcf2coanAF(full.vcf)

#Genotype data
out.name<-"coancestry_gty.txt"
gts<-extract.gt.vcf(vcf)
co.gt<-do.call(cbind,apply(gts,1,function(gt){
  snp.name<-paste(gt["#CHROM"],trimws(gt["POS"]),sep=".")
  rname<-paste(snp.name,gt["REF"],sep="_")
  aname<-paste(snp.name,gt["ALT"],sep="_")
  g<-do.call(rbind,lapply(gt[10:length(gt)],function(x) {
    strsplit(as.character(x),"/")[[1]]}))
  g[g=="0"]<-rname
  g[g=="1"]<-aname
  g[g=="."]<-0
  return(as.data.frame(g))
}))

write.table(co.gt,out.name,row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
