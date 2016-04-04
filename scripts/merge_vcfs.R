#Author: Sarah P. Flanagan
#Date: 9 February 2016
#Purpose: Merge FORMAT files exported from vcf.

setwd("E:/ubuntushare/SCA/results/biallelic")
gt<-FALSE
if(gt==TRUE){
vcf1<-read.delim("biallelic_maternal.GT.FORMAT")
vcf2<-read.delim("fem.GT.FORMAT")
vcf1$index<-paste(vcf1$CHROM,vcf1$POS,sep=".")
vcf2$index<-paste(vcf2$CHROM,vcf2$POS,sep=".")
vcf<-merge(vcf1,vcf2, by="index")

#remove the oddly duplicated ones
addedon<-vcf[duplicated(vcf$index),"index"]
if(!is.null(dim(addedon))) vcf<-vcf[!(vcf$index %in% addedon),]

#otherwise, clean it up and write to file.
drops<-c("index","CHROM.y","POS.y")
vcf<-vcf[,!(colnames(vcf) %in% drops)]
colnames(vcf)[1:2]<-c("CHROM","POS")
write.table(vcf,"biallelic.gt.vcf",col.names=T,row.names=F,
	quote=F,sep='\t')
} else{
vcf1<-read.delim("biallelic_maternal.vcf",comment.char="#",sep='\t',header=F)
vcf2<-read.delim("../stacks/batch_1.vcf",comment.char="#",sep='\t',header=F)

header.start<-grep("#CHROM",scan("biallelic_maternal.vcf",what="character"))
header1<-scan("biallelic_maternal.vcf",what="character")[header.start:
	(header.start+ncol(vcf1)-1)]
colnames(vcf1)<-header1

	
header.start<-grep("#CHROM",scan("../stacks/batch_1.vcf",what="character"))
header2<-scan("../stacks/batch_1.vcf",what="character")[header.start:
	(header.start+ncol(vcf2)-1)]
colnames(vcf2)<-header2
if(length(strsplit(as.character(vcf2[1,10]),":")[[1]])>1){
	new<-vcf2[,1:3]
	for(i in 10:ncol(vcf2)){
	new<-cbind(new,
		sapply(vcf2[,i],function(x) {
			strsplit(as.character(x),":")[[1]][1]})
		)
	}
	colnames(new)<-colnames(vcf2[,c(1:3,10:ncol(vcf2))])
	vcf2<-new
}
#remove the mom whose dad/off combo isn't right:  MOMple_PRM077_align
remove.name<-"MOMple_PRM077_align"
vcf1<-vcf1[,!(colnames(vcf1) %in% remove.name)]
#need to remove individuals in both vcf files

col.id<-c(colnames(vcf1)[1:3],colnames(vcf1)[!(colnames(vcf1) %in% 
	colnames(vcf2))])
vcf1<-vcf1[,col.id]
vcf1$index<-paste(vcf1$`#CHROM`,vcf1$ID,vcf1$POS,sep=".")
vcf2$index<-paste(vcf2$`#CHROM`,vcf1$ID,vcf2$POS,sep=".")
vcf<-merge(vcf1,vcf2, by="index")
addedon<-vcf[duplicated(vcf$index),"index"]
if(!is.null(dim(addedon))) vcf<-vcf[!(vcf$index %in% addedon),]
drops<-c("index","#CHROM.y","POS.y","ID.y")
vcf<-vcf[,!(colnames(vcf) %in% drops)]
colnames(vcf)[1:3]<-c("CHROM","POS","ID")
write.table(vcf,"biallelic_merge.vcf",col.names=T,row.names=F,
	quote=F,sep='\t')

}

