#May 7, 2018
#Sarah P. Flanagan
#Convert vcf to COLONY input

setwd("C:/Users/sflan/Documents/GitHub/SCA/results")

source("../../gwscaR/R/gwscaR.R")
source("../../gwscaR/R/gwscaR_plot.R")
source("../../gwscaR/R/gwscaR_utility.R")
source("../../gwscaR/R/gwscaR_fsts.R")
source("../../gwscaR/R/gwscaR_popgen.R")

# get the data 
vcf<-parse.vcf("drad.sub.vcf")

# calculate missingness at each locus
gts<-extract.gt.vcf(vcf)
prop.miss<-apply(gts[,10:ncol(gts)],1,function(gt){
  n<-length(gt[gt=="./."])
  prop.missing<-n/length(gt)
  return(prop.missing)
})


# Marker types
vcf<-vcf[which(prop.miss<=0.025),] #1489 loci
marker.types<-rbind(vcf$ID,rep(0,length(vcf$ID)),rep(0.2,length(vcf$ID)),rep(0.04,length(vcf$ID)))
write.table(marker.types,"colony/marker.types1489.txt",col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)

# Change the format and shorten names
co.gt<-vcf2coanGT(vcf,"colony/all_gty.txt") #this puts it in the correct format
rownames(co.gt)<-gsub("sample_(\\w+.*)_align","\\1",rownames(co.gt))

# Offspring genotypes
off.gt<-co.gt[grep("OFF",rownames(co.gt)),] #160
write.table(off.gt,"colony/off_gty.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)

# Male genotypes
mal.gt<-co.gt[c(grep("PRM",rownames(co.gt)),grep("NPM",rownames(co.gt))),] #167
write.table(mal.gt,"colony/mal_gty.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)

# Female genotypes
fem.gt<-co.gt[grep("FEM",rownames(co.gt)),] #57
write.table(fem.gt,"colony/fem_gty.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)


# Known paternity
dadIDs<-grep("PRM",rownames(co.gt),value = TRUE)
KnownPat<-lapply(dadIDs,function(id,offIDs){
  offID<-grep(gsub("PRM(\\d+)","\\1",id),offIDs,value=TRUE)
  if(id=="PRM086R"){ offID<-"OFF086" }
  if(id=="PRM086-23"){ offID<-"OFF08623" }
  known<-c(id,offID)
  if(length(known)>1){
    suppressWarnings(write.table(t(known),"colony/knownPat.txt",sep='\t',
                                 append = TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE))
  }
},offIDs=rownames(off.gt))

# Excluded paternal sibships
offIDs<-rownames(off.gt)
PatExc<-lapply(offIDs,function(id,allIDs){
  if(length(grep("-",id))>0){ nid<-gsub("(\\d+)-\\d+","\\1",id) 
  }else{ nid<-id }
  excluded<-c(id,allIDs[-which(allIDs%in%c(id,nid))])
  suppressWarnings(write.table(t(excluded),"colony/exclPatSib.txt",sep='\t',
                               append = TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE))
},allIDs=offIDs)

##### ANALYZE COLONY INITIAL #####

col.pat<-read.delim("C:/ZSL/Colony/parentage/parentage.Paternity",sep="")
col.mat<-read.delim("C:/ZSL/Colony/parentage/parentage.Maternity",sep="")

# females
avg.fem.mat<-nrow(col.mat)/length(unique(col.mat$InferredMum1[col.mat$ProbMum1>0.9]))

# males
nrow(col.pat[gsub("OFF(\\d+)(-\\d+)?","\\1",col.pat$OffspringID)==gsub("\\w{3}(\\d+)(.*)?","\\1",col.pat$InferredDad1),])/nrow(col.pat)
dup.dads<-col.pat$InferredDad1[duplicated(col.pat$InferredDad1)]
col.pat[col.pat$InferredDad1 %in% dup.dads,]

# check the ones that actually have dad in dataset
KnownPat<-read.delim("colony/knownPat.txt",header=FALSE)
known<-col.pat[col.pat$OffspringID %in% KnownPat$V2,]
nrow(known[gsub("OFF(\\d+)(-\\d+)?","\\1",known$OffspringID)==gsub("\\w{3}(\\d+)(.*)?","\\1",known$InferredDad1),])/nrow(known)

known[gsub("OFF(\\d+)(-\\d+)?","\\1",known$OffspringID)!=gsub("\\w{3}(\\d+)(.*)?","\\1",known$InferredDad1),]
nrow(known[gsub("OFF(\\d+)(-\\d+)?","\\1",known$OffspringID)==gsub("\\w{3}(\\d+)(.*)?","\\1",known$InferredDad1),])/nrow(known)

# Error rates
err<-read.delim("C:/ZSL/Colony/parentage/parentage.ErrorRate",sep=",")
col.keep.loci<-err[err$DropRateEst<0.02,]

## Re-run with only loci in col.keep.loci
vcf<-vcf[vcf$ID %in% col.keep.loci$MarkerID,]
marker.types<-rbind(vcf$ID,rep(0,length(vcf$ID)),rep(0.2,length(vcf$ID)),rep(0.04,length(vcf$ID)))
write.table(marker.types,"colony/marker.types726.txt",col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)

# Change the format and shorten names
co.gt<-vcf2coanGT(vcf,"colony/all_gty726.txt") #this puts it in the correct format
rownames(co.gt)<-gsub("sample_(\\w+.*)_align","\\1",rownames(co.gt))

# Offspring genotypes
off.gt<-co.gt[grep("OFF",rownames(co.gt)),] #160
write.table(off.gt,"colony/off_gty726.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)

# Male genotypes
mal.gt<-co.gt[c(grep("PRM",rownames(co.gt)),grep("NPM",rownames(co.gt))),] #167
write.table(mal.gt,"colony/mal_gty726.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)

# Female genotypes
fem.gt<-co.gt[grep("FEM",rownames(co.gt)),] #57
write.table(fem.gt,"colony/fem_gty726.txt",row.names = TRUE,col.names = FALSE,sep='\t',quote=FALSE)


