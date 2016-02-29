#Author: Sarah P. Flanagan
#Date: 26 February 2016
#Purpose: Prune loci to include only those RAD loci found in 75% of the
#original RAD seq dataset and in 75% of the ddRAD seq dataset

setwd("../results/stacks")
vcf<-read.table("batch_1.vcf")
#get the header
colnames(vcf)<-read.table("batch_1.vcf.header.txt",comment.char="",header=F,
	stringsAsFactors=F)
info.fields<-colnames(vcf)[seq(1,9)]
dd.ind<-read.table("../biallelic/ind.info.ddRAD.txt",header=F)
or.ind<-read.table("../biallelic/ind.info.psti.txt",header=F)

dd.vcf<-vcf[,colnames(vcf) %in% c(info.fields,as.character(dd.ind$V1))]
or.vcf<-vcf[,colnames(vcf) %in% c(info.fields,as.character(or.ind$V1))]

dd.prop<-apply(dd.vcf,1,function(x){
	total<-length(x)-9
	empty<-length(unlist(lapply(x,grep,"./.:0:.,.:.,.,.")))
	prop<-empty/total
	return(prop) })
ddk.vcf<-dd.vcf[dd.prop<0.25,]
ddk.vcf$Locus<-paste(ddk.vcf[,1],".",ddk.vcf[,2],sep="")


or.prop<-apply(or.vcf,1,function(x){
	total<-length(x)-9
	empty<-length(unlist(lapply(x,grep,"./.:0:.,.:.,.,.")))
	prop<-empty/total
	return(prop) })
ork.vcf<-or.vcf[or.prop<0.25,]
ork.vcf$Locus<-paste(ork.vcf[,1],".",ork.vcf[,2],sep="")

keep.locus<-ddk.vcf[ddk.vcf$Locus %in% ork.vcf$Locus,]
write.table(keep.locus[,c("#CHROM","POS","ID","Locus")],
	"LociToKeep.txt",sep="\t",row.names=F,quote=F,
	col.names=c("Scaffold","Position","CatalogID","LocusID"))
write.table(keep.locus[,"ID"],
	"CatalogIDsToKeep.txt",sep="\t",row.names=F,quote=F,
	col.names=F)

#prune the original vcf
vcf$Locus<-paste(vcf[,1],".",vcf[,2],sep="")
vcf.prune<-vcf[vcf$Locus %in% keep.locus$Locus,]
write.table(subset(vcf.prune,select=-Locus),"batch_1.pruned.vcf",sep="\t",
	col.names=T,row.names=F,quote=F)


