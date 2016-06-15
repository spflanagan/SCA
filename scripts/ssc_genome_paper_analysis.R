#Author: Sarah P. Flanagan
#Last Updated: 27 April 2016
#Start Date: 27 April 2016
#Purpose: Conduct a linkage disequilibrium analysis

rm(list=ls())
setwd("E:/ubuntushare/SCA/results/genome_paper")
library(RColorBrewer)
library(gplots)
library(qvalue)
source("E:/ubuntushare/SCA/scripts/plotting_functions.R")

lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
##############################################################################
#####Fsts
##############################################################################
fm.plot<-read.delim("../biallelic/fm.plot.txt")
fm.plot<-fm.plot[,c("Chrom","Pos","LocID","FEM.MAL","Locus")]
fm.plot<-fm.plot[order(fm.plot$FEM.MAL),]#ascending
fm.top1<-fm.plot[round(nrow(fm.plot)*0.99),"FEM.MAL"]
fm.out1<-fm.plot[fm.plot$FEM.MAL >= fm.top1,]


png("male-female.png",height=100,width=300,units="mm",res=300)
par(oma=c(2,2,2,2),mar=c(4,0,0,0))
fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0,bp.name="Pos",sig.col=c("green4","black"))
axis(2,at=c(0,0.1,0.2,0.3),pos=0,las=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0,line=1)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fm[fm$Chrom ==lgs[i],"Pos"]),y=-0.002,
		labels=lgs[i], srt=45, adj=1, xpd=TRUE)
	last<-max(fm[fm$Chrom ==lgs[i],"Pos"])
}
dev.off()

#####STACKS FSTS
stacks<-read.delim("../stacks/batch_1.fst_FEM-PRM.tsv",header=T)
outpoints<-stacks[order(stacks$Fst),]
stacks.top1<-outpoints[round(nrow(outpoints)*0.99),"Fst"]
#outpoints<-stacks[order(stacks$Corrected.Fst),]
#corr.top1<-outpoints[round(nrow(outpoints)*0.99),"Fst"]
#different Fst versions don't make a difference

png("male-female_stacks.png",height=100,width=300,units="mm",res=300)
par(oma=c(2,2,2,2),mar=c(4,0,0,0))
sp<-fst.plot(stacks, ci.dat=c(stacks.top1,0),fst.name="Fst", 
	chrom.name="Chr", axis.size=0,bp.name="BP",sig.col=c("green4","grey7"))
axis(2,at=seq(-0.2,0.8,0.2),pos=0,las=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0,line=1)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(sp[sp$Chr ==lgs[i],"BP"]),y=-0.202,
		labels=lgs[i], srt=45, adj=1, xpd=TRUE)
	last<-max(sp[sp$Chr ==lgs[i],"BP"])
}
dev.off()


##############################################################################
#####SNPSTATS
##############################################################################
snpstats<-read.table("snpstats1_out.txt",header=T)
snp.plots<-snpstats[,c("LG","pos","p_val.3")]
snp.plots$plotp<-(log10(snp.plots$p_val.3))*-1
snp.plots<-snp.plots[snp.plots$plotp != "Inf" &snp.plots$plotp != "-Inf",]
lfdr(snp.plots$p_val.3)#doesn't work, p-values out of range.

png("genome_rad_fig1.png", width=10,height=7,units="in",res=300)
par(mfrow=c(2,1),oma=c(2,2,2,2),mar=c(1,0,1,0))
fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0,bp.name="Pos",sig.col=c("green4","black"))
axis(2,at=c(0,0.1,0.2,0.3),pos=0,las=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=F, cex=1,las=0,line=1)
text(x=1900,y=0.205,"A")

par(mar=c(3,0,0,0))
g<-fst.plot(snp.plots,ci.dat=c(0,0),sig.col=c("black","black"),
	fst.name="plotp",chrom.name="LG",bp.name="pos",axis.size=0)
axis(2,at=seq(-1.6,4.6,0.4),pos=0,las=1)
mtext(expression(-log[10]*P), 2, outer=F, cex=1,las=0,line=1)
text(x=2200,y=4.2,"B")

last<-0
for(i in 1:length(lgs)){
	text(x=mean(g[g$LG ==lgs[i],"pos"]),y=-1.8,
		labels=lgs[i], srt=45, adj=1, xpd=TRUE)
	last<-max(g[g$LG ==lgs[i],"pos"])
}
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
ld.files<-list.files(pattern="adults_maf1_LG\\d+.txt")
ld.files<-c(ld.files,"ld_matrix_adults_maf1_lg1.txt")
ld.order<-c("lg1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10",
	"LG11","LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20",
	"LG21","LG22")
for(i in 1:length(ld.order)){
	filename<-ld.files[gsub("ld_matrix_adults_maf1_([A-z]+\\d+).txt","\\1",
		ld.files) %in% ld.order[i]]
	ld.dat<-read.table(filename,header=T,row.names=1)
	ld.heatmap(ld.dat,
		name=gsub("(ld_matrix_adults_maf1_[A-z]+\\d+).txt","\\1.png",filename))
}

###FINAL FIGURE PLOTTED USING GIMP2
ld.dat<-read.table("ld_matrix_adults_maf1_LG12.txt",header=T,row.names=1)
sumstats<-read.delim("../stacks/batch_1.sumstats.tsv",sep='\t',header=T,skip=4)
sumstats$locus<-paste(sumstats$Chr,sumstats$BP,sumstats$Locus.ID,sep=".")
lg12.sum<-sumstats[sumstats$Chr %in% "LG12",]
lg12.ld<-sumstats[sumstats$locus %in% rownames(ld.dat),]
lg12.prm<-lg12.ld[lg12.ld$Pop.ID == "PRM",]
lg12.fem<-lg12.ld[lg12.ld$Pop.ID == "FEM",]
large.ld<-rownames(ld.dat)[rowSums(ld.dat) > 1000]

large.prm.sum<-lg12.prm[lg12.prm$locus %in% large.ld,]
others.prm.sum<-lg12.prm[!(lg12.prm$locus %in% large.ld),]
large.fem.sum<-lg12.fem[lg12.fem$locus %in% large.ld,]
others.fem.sum<-lg12.fem[!(lg12.fem$locus %in% large.ld),]

png("elevatedLD_sumstats.png",height=10,width=7,units="in",res=300)
par(mfrow=c(4,2),oma=c(2,2,2,2),mar=c(2,2,2,2))
hist(others.prm.sum$P,col=alpha("grey",0.5),breaks=seq(0,1,0.1),
	main="Non-Elevated Loci",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(others.fem.sum$P,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.1),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Major Allele Frequency",1,line=1.5,cex=0.75)
legend("topleft",c("Males","Females"),pch=15,bty='n',
	col=c(alpha("grey",0.5),alpha("dodgerblue",0.5)))
hist(large.prm.sum$P,col=alpha("grey",0.5),breaks=seq(0,1,0.1),
	main="Elevated LD Loci",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(large.fem.sum$P,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.1),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Major Allele Frequency",1,line=1.5,cex=0.75)

hist(others.prm.sum$Pi,col=alpha("grey",0.5),breaks=seq(0,1,0.1),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(others.fem.sum$Pi,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.1),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Nucleotide Diversity",1,line=1.5,cex=0.75)
hist(large.prm.sum$Pi,col=alpha("grey",0.5),breaks=seq(0,1,0.1),ylim=c(0,60),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(large.fem.sum$Pi,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.1),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Nucleotide Diversity",1,line=1.5,cex=0.75)

hist(others.prm.sum$Fis,col=alpha("grey",0.5),breaks=seq(-1,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=-1,las=1)
hist(others.fem.sum$Fis,col=alpha("dodgerblue",0.5),add=T,breaks=seq(-1,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=-1,las=1)
mtext(expression(italic(F)[IS]),1,line=1.5,cex=0.75)
hist(large.prm.sum$Fis,col=alpha("grey",0.5),breaks=seq(-1,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=-1,las=1)
hist(large.fem.sum$Fis,col=alpha("dodgerblue",0.5),add=T,breaks=seq(-1,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=-1,las=1)
mtext(expression(italic(F)[IS]),1,line=1.5,cex=0.75)

hist(others.prm.sum$Obs.Het,col=alpha("grey",0.5),breaks=seq(0,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(others.fem.sum$Obs.Het,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Observed Heterozygosity",1,line=1.5,cex=0.75)
hist(large.prm.sum$Obs.Het,col=alpha("grey",0.5),breaks=seq(0,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
hist(large.fem.sum$Obs.Het,col=alpha("dodgerblue",0.5),add=T,breaks=seq(0,1,0.2),
	main="",xaxt='n',yaxt='n',xlab="",ylab="")
axis(1,pos=0)
axis(2,pos=0,las=1)
mtext("Observed Heterozygosity",1,line=1.5,cex=0.75)

dev.off()

ld.heatmap(ld.dat,make.file=F)
ld.large<-ld.dat[rownames(ld.dat) %in% large.ld,colnames(ld.dat) %in% large.ld]
ld.heatmap(ld.large,make.file=F)

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


