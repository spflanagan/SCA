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
#indnames2<-colnames(vcf2)[4:ncol(vcf2)]
#indnames2<-gsub("sample_(\\w+\\d+.*)_align","\\1",indnames2)
#colnames(vcf2)[4:ncol(vcf2)]<-indnames2
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
ind.info<-data.frame(filenameID=colnames(vcf)[10:ncol(vcf)],
	ID=colnames(vcf)[10:ncol(vcf)],
	sex=substr(colnames(vcf)[10:ncol(vcf)],1,3),
	age=substr(colnames(vcf)[10:ncol(vcf)],1,3),
	status="0",
	phen=substr(colnames(vcf)[10:ncol(vcf)],1,3),stringsAsFactors=F)

#set the sexes
ind.info$sex[ind.info$sex=="MOM"]<-"0"
ind.info$sex[ind.info$sex=="NPM" | ind.info$sex=="PRM"]<-"MAL"
ind.info$sex[ind.info$sex=="OFF"]<-"0"

#set age
ind.info$age[ind.info$age=="MOM"]<-"0"
ind.info$age[ind.info$age=="NPM" | ind.info$age=="PRM" | 
	ind.info$age=="FEM"]<-"ADULT"
ind.info$age[ind.info$age=="OFF"]<-"JUVIE"

#set phenotype
ind.info$phen[ind.info$phen=="NPM" | ind.info$phen=="PRM" | 
	ind.info$phen=="FEM" | ind.info$phen=="OFF"]<-"0"
#if we wanted to include pregnant/nonpregnant, this is where we would set that

#write.table(ind.info,"ind_info_vcf.txt",quote=F,row.names=F,col.names=F,
#	sep='\t')
	

}

