#Author: Sarah P. Flanagan
#Date: 28 March 2016
#Purpose: Parse and prune blast output

prune.blast<-function(blast.table){
	loc.bn.split<-split(blast.table,blast.table$qseqid)
	loc.bn.keep<-NULL
	for(i in 1:length(loc.bn.split)){
	#by length
	match.length<-loc.bn.split[[i]]$send-loc.bn.split[[i]]$sstart
	ok<-loc.bn.split[[i]][abs(match.length)>50,]
	#by percent identity
	keep<-which.max(ok$pident)
	if(length(keep)>0) {
		loc.bn.keep<-rbind(loc.bn.keep,loc.bn.split[[i]][keep,]) }
	}
	return(loc.bn.keep)
}

setwd("E:/ubuntushare/SCA/results/biallelic_outliers/rad_locus")
locus.blast<-read.table("top1.out.radloc.blastn",sep='\t')
colnames(locus.blast)<-c("qseqid","sseqid","pident","length","mismatch",
	"gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle")
loc.bn<-prune.blast(locus.blast)
write.table(loc.bn,"top1.radloc.keep.blastn",sep='\t',quote=F,row.names=F)


setwd("E:/ubuntushare/SCA/results/biallelic_outliers/rad_region")
region.blast<-read.table("merged_extracts.blastn",sep='\t')
colnames(region.blast)<-c("qseqid","sseqid","pident","length","mismatch",
	"gapopen","qstart","qend","sstart","send","evalue","bitscore","stite")
reg.bn<-prune.blast(region.blast)
write.table(reg.bn,"merged_extracts.keep.blstn",sep='\t',quote=F,row.names=F)