#Author: Sarah P. Flanagan
#Date: 24 March 2016
#Purpose: Analyze the alignments

setwd("E:/ubuntushare/SCA/results/biallelic_outliers/align_to_annotated")
sam.names<-c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT",
	"TLEN","SEQ","QUAL","OPT1","OPT2","OPT3","OPT4","OPT5","OPT6","OPT7",
	"OPT8","OPT9") 
	#POS is leftmost mapping position (1-based)
	#RNAME is reference name
gff<-read.delim("E:/ubuntushare/scovelli_genome/annotated_scovelli_scaffolds_March2016.gff",
	header=F)
gff.names<-c("sequence","source","feature","start","end","score","strand",
	"frame","attributes")
colnames(gff)<-gff.names
	#sequence is name of sequence
	#start is first bp(one base offset)
	#end is last bp (one base offset)

sam2gff<-function(samrow){
	rel.gff<-gff[gff$sequence %in% samrow["RNAME"],]
	out<-rel.gff[as.numeric(rel.gff$start) <= as.numeric(samrow["POS"]) & 
		as.numeric(rel.gff$end) >= as.numeric(samrow["POS"]) &
		 rel.gff$feature != "contig",]
	if(nrow(out)>0){
		out<-cbind(QNAME=rep(samrow["QNAME"],nrow(out)), 
			POS=rep(samrow["POS"],nrow(out)), out)
		return(out)
	}
}

#sam<-read.delim("top1.out.radloc_align.sam",comment.char="@",header=F)
#colnames(sam)<-sam.names
#matches<-do.call("rbind",apply(sam,1,sam2gff))
#write.table(matches,"sam_to_gff.txt",sep="\t",quote=F,col.names=T,row.names=F)

#aj.sam<-read.delim("aj_outlier_loci_align.sam",comment.char="@",header=F)
#colnames(aj.sam)<-sam.names[1:ncol(aj.sam)]
#aj.matches<-do.call("rbind",apply(aj.sam,1,sam2gff))
#write.table(aj.matches,"aj_sam_to_gff.txt",sep="\t",quote=F,col.names=T,row.names=F)

fm.sam<-read.delim("fm_1outlier_loci_align.sam",comment.char="@",header=F)
colnames(fm.sam)<-sam.names[1:ncol(fm.sam)]
fm.matches<-do.call("rbind",apply(fm.sam,1,sam2gff))
write.table(fm.matches,"fm_1sam_to_gff.txt",sep="\t",quote=F,col.names=T,row.names=F)

mo.sam<-read.delim("mo_1outlier_loci_align.sam",comment.char="@",header=F)
colnames(mo.sam)<-sam.names[1:ncol(mo.sam)]
mo.matches<-do.call("rbind",apply(mo.sam,1,sam2gff))
write.table(mo.matches,"mo_1sam_to_gff.txt",sep="\t",quote=F,col.names=T,row.names=F)

pj.sam<-read.delim("pj_1outlier_loci_align.sam",comment.char="@",header=F)
colnames(pj.sam)<-sam.names[1:ncol(pj.sam)]
pj.matches<-do.call("rbind",apply(pj.sam,1,sam2gff))
write.table(pj.matches,"pj_1sam_to_gff.txt",sep="\t",quote=F,col.names=T,row.names=F)

