#Author: Sarah P. Flanagan
#Last Updated: 27 April 2016
#Start Date: 27 April 2016
#Purpose: Conduct a linkage disequilibrium analysis

rm(list=ls())
setwd("E:/ubuntushare/SCA/results/ld")
library(RColorBrewer)
library(gplots)

##############################################################################
#####SNPSTATS
##############################################################################
snpstats<-read.table("../sexlinked/snpstats1_out.txt",header=T)
snp.plots<-snpstats[,c("LG","pos","p_val.3")]
snp.plots$plotp<--log10(snp.plots$p_val.3)
snp.plots<-snp.plots[snp.plots$plotp != "Inf" &snp.plots$plotp != "-Inf",]
png("logP_G.png", width=10,height=7,units="in",res=300)
g<-fst.plot(snp.plots,ci.dat=c(0,0),sig.col=c("black","black"),
	fst.name="plotp",chrom.name="LG",bp.name="pos")
dev.off()


##############################################################################
#####LD ANALYSIS
##############################################################################
ld.heatmap<-function(ld.matrix, name="ld.heatmap",pdf=F,make.file=T){
	if(nrow(ld.matrix)%%20 != 0){
		nmax<-nrow(ld.matrix)-(nrow(ld.matrix)%%20)
		starts<-c(seq(1,nmax+1,20))
		ends<-c(seq(20,nmax,20),nrow(ld.matrix))
	}
	else {
		starts<-starts<-c(seq(1,nrow(ld.matrix)-19,20))
		ends<-c(seq(20,nrow(ld.matrix),20))

	}
	heatcolors<-colorRampPalette(c("white","yellow","red"))(n=200)

	if(make.file==T){
		if(pdf==T) { pdf(paste(name,"pdf",sep="."),height=10,width=10) } 
		else {
			png(paste(name,"png",sep="."),height=10,width=10,units="in",
				res=300) }
	}
	par(mfrow=c(length(starts),length(starts)),oma=c(0,0,0,0),mar=c(0,0,0,0))
	for(i in 1:length(starts)){
		for(j in 1:length(ends)){
			#heatmap.2(as.matrix(ld.matrix[
			#	starts[i]:ends[i],starts[j]:ends[j]]),
			#	dendrogram="none",tracecol="NA",labCol="",labRow="",key=F,
			#	Colv=F,Rowv=F,lwid=c(0.5,4),lhei=c(0.5,4),new=F,
			#	col=heatcolors)
			image(as.matrix(ld.dat[
					starts[i]:ends[i],starts[j]:ends[j]]),
				xaxt='n',yaxt='n',col=heatcolors,bty='n')
			#if(i==j) print(paste(starts[i],":",ends[i],",",starts[j],":",
			#	ends[j],sep=""))
	}}
	if(make.file==T) dev.off()
}
ld.files<-list.files(pattern="LG\\d+.txt")
ld.files<-c(ld.files,"ld_matrix_lg1.txt")
ld.order<-c("lg1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10",
	"LG11","LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20",
	"LG21","LG22")
for(i in 1:length(ld.order)){
	filename<-ld.files[gsub("ld_matrix_([A-z]+\\d+).txt","\\1",
		ld.files) %in% ld.order[i]]
	ld.dat<-read.table(filename,header=T,row.names=1)
	ld.heatmap(ld.dat,
		name=gsub("(ld_matrix_[A-z]+\\d+).txt","\\1.png",filename))
}


#############################################################################
#vcf<-read.delim("../stacks/batch_1.vcf",comment.char="#",sep='\t')
#header<-scan("../stacks/batch_1.vcf",what="character")[header.start:
#	(header.start+ncol(vcf)-1)]
#colnames(vcf)<-header1
#if(length(strsplit(as.character(vcf[1,10]),":")[[1]])>1){
#	new<-vcf[,1:3]
#	for(i in 10:ncol(vcf)){
#	new<-cbind(new,
#		sapply(vcf[,i],function(x) {
#			strsplit(as.character(x),":")[[1]][1]})
#		)
#	}
#	colnames(new)<-colnames(vcf[,c(1:3,10:ncol(vcf))])
#	vcf<-new
#}
#colnames(vcf)<-header[c(1:3,10:length(header))]
#vcf.chrom<-split(vcf,vcf$`#CHROM`)

#calc.ld<-function(vcf.list,row1,row2){
#	loc1<-vcf.list[row1,4:ncol(vcf.list)]
#	loc2<-vcf.list[row2,4:ncol(vcf.list)]
#	joint<-rbind(loc1,loc2)
#	joint<-joint[,joint[1,] != "./." & joint[2,]!="./."]
#	joint.freqs<-data.frame("0"=c(0,0),"1"=c(0,0),row.names=c("0","1"))
#	freqs1<-data.frame("0"=0,"1"=0)
#	freqs2<-data.frame("0"=0,"1"=0)
#	for(i in 1:ncol(joint)){
#		mat1<-as.numeric(strsplit(as.character(joint[1,i]),"/")[[1]][1])+1
#		pat1<-as.numeric(strsplit(as.character(joint[1,i]),"/")[[1]][2])+1
#		mat2<-as.numeric(strsplit(as.character(joint[2,i]),"/")[[1]][1])+1
#		pat2<-as.numeric(strsplit(as.character(joint[2,i]),"/")[[1]][2])+1
#		freqs1[mat1]<-freqs1[mat1]+1
#		freqs1[pat1]<-freqs1[pat1]+1
#		freqs2[mat2]<-freqs2[mat2]+1
#		freqs2[pat2]<-freqs2[pat2]+1
#		joint.freqs[mat1,mat2]<-joint.freqs[mat1,mat2]+1
#		joint.freqs[pat1,pat2]<-joint.freqs[pat1,pat2]+1
#	}
#	joint.freqs<-joint.freqs/(2*ncol(joint))
#	freqs1<-freqs1/(2*ncol(joint))
#	freqs2<-freqs2/(2*ncol(joint))
#	
#	d<-data.frame("0"=c(0,0),"1"=c(0,0),row.names=c("0","1"))
#	dmax<-d
#	for(i in 1:2){
#		for(j in 1:2){	
#			d[i,j]<-joint.freqs[i,j]-(freqs1[i]*freqs2[j])
#			if(d[i,j]<0){	
#				dm<-min((freqs1[i]*freqs2[j]),
#					((1-freqs1[i])*(1-freqs2[j])))
#			}else{
#				dm<-min(((1-freqs1[i])*freqs2[j]),
#					((freqs1[i])*(1-freqs2[j])))
#			}
#			dmax[i,j]<-dm
#		}
#	}
#	dprime<-0
#	for(i in 1:2){
#		for(j in 1:2){
#			if(freqs1[i] > 0 & freqs2[j] > 0){
#				if(dmax[i,j]>0){
#					dprime<-dprime+(freqs[i]*freqs[j]*
#						abs(d[i,j])/dmax[i,j])
#				} else { dprime<- -5 }
#			}
#		}
#	}
#	return(dprime)
#}
#ld<-list(rep(data.frame(),length(vcf.chrom)))
#for(f in 1:length(vcf.chrom)){
#	for(ff in 1:length(vcf.chrom[[f]])){
#		for(fff in 1:length(vcf.chrom[[f]])){
#			ld[[f]][ff,fff]<-calc.ld(vcf.chrom[[f]],ff,fff)
#		}
#	}
#
#}


