#Author: Sarah P. Flanagan
#Last Updated: 27 April 2016
#Started Date: 11 February 2016
#Purpose: Analyze the output from gwsca_biallelic.

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")
lgs<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7","LG8","LG9","LG10","LG11",
	"LG12","LG13","LG14","LG15","LG16","LG17","LG18","LG19","LG20","LG21",
	"LG22")
lgn<-seq(1,22)

#####################################FILES###################################
gw.fst<-read.delim("gwsca_fsts_new.txt")
	gw.fst$Locus<-paste(gw.fst$Chrom,gw.fst$LocID,gw.fst$Pos,sep=".")
gw.sum<-read.delim("gwsca_summary_new.txt")
	gw.sum$Locus<-paste(gw.sum$Chrom,gw.sum$LocID,gw.sum$Pos,sep=".")
aj.prune<-read.table("aj.plot.txt",header=T,sep='\t')
fm.prune<-read.table("fm.plot.txt",header=T,sep='\t')
mo.prune<-read.table("mo.plot.txt",header=T,sep='\t')
hd<-read.table(sep="\t", header=F,
	"../monnahan/ML.2016.pipefish_input_full.txt")
	colnames(hd)<-c("chrom","pos","parent","q_no","mu_no","sig_no","LL1",
	"qM_0","qF_0","mu_0","sig_0","LL0","LRT0","q_fec","mu1_fec",
	"mu2_fec","mu3_fec","sig_fec","LL2","LRT2","q_sex","qX_sex","mu_sex",
	"sig_sex","LL3","LRT3")
snps<-read.table("../stacks/batch_1.catalog.snps.tsv",sep='\t',comment.char="#")
	colnames(snps)<-c("SqlID","SampleID","LocusID","Column","Type","LR","Rank1",
	"Rank2","Rank3","Rank4")
tags<-read.table("../stacks/batch_1.catalog.tags.tsv",sep='\t',comment.char="#")
	colnames(tags)<-c("SqlID","SampleID","LocusID","Chrom","BP","Strand",
	"SeqType","StackComponent","SeqID","Seq","Deleveraged","Blacklisted",
	"Lumberjackstack","LogLike")
vcf<-read.delim("../stacks/batch_1.vcf",comment.char="#",sep='\t',header=F)
	header.start<-grep("#CHROM",scan("../stacks/batch_1.vcf",what="character"))
	header<-scan("../stacks/batch_1.vcf",what="character")[header.start:
	(header.start+ncol(vcf)-1)]
	colnames(vcf)<-header
#############################################################################

#################################FUNCTIONS###################################
extract.info<-function(out.sum,tags.dat,snps.dat){
	ex.tag<-tags[tags$LocusID %in% out.sum$LocID,c("LocusID","Chrom","BP","Seq")]
	ex.snp<-snps[snps$LocusID %in% out.sum$LocID,]
	ex.snp.split<-split(ex.snp,ex.snp$LocusID)
	ex.dat<-do.call(rbind,lapply(ex.snp.split,function(x){
		num<-nrow(x)
		ranks<-list(as.character(x$Rank1[1]),as.character(x$Rank2[1]),
			as.character(x$Rank3[1]),as.character(x$Rank4[1]))
		for(i in 2:num){
			ranks[1]<-paste(ranks[1],x$Rank1[i],sep=",")
			ranks[2]<-paste(ranks[2],x$Rank2[i],sep=",")
			ranks[3]<-paste(ranks[3],x$Rank3[i],sep=",")
			ranks[4]<-paste(ranks[4],x$Rank4[i],sep=",")
		}
		df<-data.frame(LocusID=levels(factor(x$LocusID)),NumSNPs=num,
			do.call(cbind,ranks))
		colnames(df)<-c("LocusID","NumSNPs","SNP1","SNP2","SNP3","SNP4")
		return(df)
	}))
	ex.dat<-merge(ex.dat,ex.tag,by="LocusID")
	ex.split<-split(out.sum,out.sum$LocID)
	ex.sum.m<-do.call(rbind,lapply(ex.split,function(x){
		y<-split(x,factor(x$Pop))
		num<-nrow(y[[1]])
		if(num>1){
			yinfo<-lapply(y,function(z){
				info<-list(z$Allele1Freq[1],z$Ho[1],z$Hs[1],z$N[1])
				for(i in 2:num){
					info[1]<-paste(info[1],z$Allele1Freq[i],sep=",")
					info[2]<-paste(info[2],z$Ho[i],sep=",")
					info[3]<-paste(info[3],z$Hs[i],sep=",")
					info[4]<-paste(info[4],z$N[i],sep=",")
				}
				zout<-do.call(cbind,info)
				colnames(zout)<-c("AF","Ho","Hs","N")
				return(zout)
			})
			posinfo<-y[[1]]$Pos[1]
			for(i in 2:num){
				posinfo<-y[[1]]$Pos[i]
			}
		} else {
			yinfo<-lapply(y,function(z){
				info<-data.frame(AF=as.character(z$Allele1Freq[1]),
					Ho=as.character(z$Ho[1]),Hs=as.character(z$Hs[1]),
					N=as.character(z$N[1]))
				return(info)
			})
			posinfo<-y[[1]]$Pos[1]
		}
		df<-data.frame(LocusID=as.factor(levels(factor(x$LocID))),
			NumSNPs=as.numeric(num),BPs=posinfo,
			do.call(cbind,yinfo))
		colnames(df)<-c("LocusID","NumSNPsOut","BPs","Pop1AF","Pop1Ho","Pop1Hs",
			"Pop1N","Pop2AF","Pop2Ho","Pop2Hs","Pop2N")
		return(df)
	}))

	ex.dat<-merge(ex.dat,ex.sum.m,by="LocusID")
	return(ex.dat)
}

parse.snps<-function(snp){
	snp.split<-split(snp,snp$LocusID)
	ex.dat<-do.call(rbind,lapply(snp.split,function(x){
		num<-nrow(x)
		ranks<-list(as.character(x$Rank1[1]),as.character(x$Rank2[1]),
			as.character(x$Rank3[1]),as.character(x$Rank4[1]))
		for(i in 2:num){
			ranks[1]<-paste(ranks[1],x$Rank1[i],sep=",")
			ranks[2]<-paste(ranks[2],x$Rank2[i],sep=",")
			ranks[3]<-paste(ranks[3],x$Rank3[i],sep=",")
			ranks[4]<-paste(ranks[4],x$Rank4[i],sep=",")
		}
		df<-data.frame(LocusID=levels(factor(x$LocusID)),NumSNPs=num,
			do.call(cbind,ranks))
		colnames(df)<-c("LocusID","NumSNPs","SNP1","SNP2","SNP3","SNP4")
		return(df)
	}))
	return(ex.dat)
}
#############################################################################

###################################PRUNING###################################
#remove any that are not polymorphic
gw.sum<-gw.sum[!is.na(gw.sum$AA),]
#test for HWE
gw.sum$AAexp<-gw.sum$Allele1Freq*gw.sum$Allele1Freq
gw.sum$aaexp<-gw.sum$Allele2Freq*gw.sum$Allele2Freq
gw.sum$Aaexp<-1-gw.sum$aaexp-gw.sum$AAexp

gw.sum$chi<-(((gw.sum$AA-gw.sum$AAexp)^2)/gw.sum$AAexp)+
	(((gw.sum$Aa-gw.sum$Aaexp)^2)/gw.sum$Aaexp)+
	(((gw.sum$aa-gw.sum$aaexp)^2)/gw.sum$aaexp)
gw.sum$chi.result<-1-pchisq(gw.sum$chi,1) #biallelic, df=1
gw.hwe<-gw.sum[gw.sum$chi.result > 0.05,]

#prune to keep only those found in most pops
sum.prune<-gw.hwe[gw.hwe$Pop=="ADULT" | gw.hwe$Pop=="JUVIE",]
sum.sum<-tapply(sum.prune$N,sum.prune$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 393]
sum.prune<-gw.hwe[gw.hwe$Locus %in% names(sum.sum),]

##prune based on representation in the groups (in 50% of individuals)
sum.list<-split(sum.prune, sum.prune$Pop)
adt.n<-sum.list$ADULT[sum.list$ADULT$N > 216 & !is.na(sum.list$ADULT$Hs),]
juv.n<-sum.list$JUVIE[sum.list$JUVIE$N > 157 & !is.na(sum.list$JUVIE$Hs),]
fem.n<-sum.list$FEM[sum.list$FEM$N > 57& !is.na(sum.list$FEM$Hs),]
mal.n<-sum.list$MAL[sum.list$MAL$N>159& !is.na(sum.list$MAL$Hs),]
mom.n<-sum.list$MOM[sum.list$MOM$N>133& !is.na(sum.list$MOM$Hs),]
prg.n<-sum.list$PREGNANT[sum.list$PREGNANT$N>150& !is.na(sum.list$PREGNANT$Hs),]

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]
prg.n<-prg.n[prg.n$Allele1Freq > 0.05 & prg.n$Allele1Freq < 0.95,]

#write these out to read in for error analysis
write.table(prg.n[,2:4],"LociInPregMales.txt",col.names=T,row.names=F,quote=F)
write.table(juv.n[,2:4],"LociInOffspring.txt",col.names=T,row.names=F,quote=F)
#comparisons
#viability
aj.prune<-gw.fst[gw.fst$Locus %in% adt.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
aj.prune<-aj.prune[aj.prune$ADULT.JUVIE>0,]
fm.prune<-gw.fst[gw.fst$Locus %in% mal.n$Locus&gw.fst$Locus %in% fem.n$Locus, ]
fm.prune<-fm.prune[fm.prune$FEM.MAL>0,]
#sexual
mo.prune<-gw.fst[gw.fst$Locus %in% fem.n$Locus&gw.fst$Locus %in% mom.n$Locus, ]
mo.prune<-mo.prune[mo.prune$MOM.FEM>0,]

#write.table(fm.prune, "fm.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
#write.table(mo.prune, "mo.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
#write.table(aj.prune,"aj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
#############################################################################

#################################ANALYSIS####################################
aj.plot<-aj.prune[order(aj.prune$ADULT.JUVIE),] #ascending
aj.top1<-aj.plot[round(nrow(aj.plot)*0.99),"ADULT.JUVIE"]
aj.out1<-aj.plot[aj.plot$ADULT.JUVIE >= aj.top1,]
fm.plot<-fm.prune[order(fm.prune$FEM.MAL),]#ascending
fm.top1<-fm.plot[round(nrow(fm.plot)*0.99),"FEM.MAL"]
fm.out1<-fm.prune[fm.prune$FEM.MAL >= fm.top1,]
mo.plot<-mo.prune[order(mo.prune$MOM.FEM),]#ascending
mo.top1<-mo.plot[round(nrow(mo.plot)*0.99),"MOM.FEM"]
mo.out1<-mo.prune[mo.prune$MOM.FEM >= mo.top1,]


#plot with the top1%
#NOT THE ACTUAL FIGURE
mo.plot$CompLoc<-paste(mo.plot$Chrom,mo.plot$Pos,sep=".")
mo<-fst.plot(mo.plot, ci.dat=c(mo.top1,0),fst.name="MOM.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"),groups=lgs)
fm.plot$CompLoc<-paste(fm.plot$Chrom,fm.plot$Pos,sep=".")
fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
aj.plot$CompLoc<-paste(aj.plot$Chrom,aj.plot$Pos,sep=".")
aj<-fst.plot(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))

mo.plot<-mo.plot[mo.plot$Chrom %in% lgs,]
mo.plot$Chrom<-factor(mo.plot$Chrom)
fm.plot<-fm.plot[fm.plot$Chrom %in% lgs,]
fm.plot$Chrom<-factor(fm.plot$Chrom)
aj.plot<-aj.plot[aj.plot$Chrom %in% lgs,]
aj.plot$Chrom<-factor(aj.plot$Chrom)

###COMPARISONS
aj.out<-aj.plot[aj.plot$ADULT.JUVIE >= aj.top1[1],]
fm.out<-fm.plot[fm.plot$FEM.MAL >= fm.top1[1],]
mo.out<-mo.plot[mo.plot$MOM.FEM >= mo.top1[1],]
aj.out.plot<-aj[aj$ADULT.JUVIE >= aj.top1[1],]
fm.out.plot<-fm[fm$FEM.MAL >= fm.top1[1],]
mo.out.plot<-mo[mo$MOM.FEM >= mo.top1[1],]

aj.unique<-aj.out[!(aj.out$Locus %in% fm.out$Locus) & 
	!(aj.out$Locus %in% mo.out$Locus),]
fm.unique<-fm.out[!(fm.out$Locus %in% aj.out$Locus) & 
	!(fm.out$Locus %in% mo.out$Locus),]
mo.unique<-mo.out[!(mo.out$Locus %in% aj.out$Locus) &
	!(mo.out$Locus %in% fm.out$Locus),]
shared.out<-aj.out[(aj.out$LocID %in% mo.out$LocID) & 
	(aj.out$LocID %in% fm.out$LocID),]

aj.fm<-aj.out[(aj.out$LocID %in% fm.out$LocID),]
aj.mo<-aj.out[(aj.out$LocID %in% mo.out$LocID),]
fm.mo<-fm.out[fm.out$LocID %in% mo.out$LocID,]
length(levels(factor(c(as.character(aj.fm$Locus),as.character(aj.mo$Locus),
	as.character(fm.mo$Locus)))))

#####COMPARE TO MONNAHAN/KELLY
hd$Locus<-paste(hd$chrom,hd$pos,sep=".")

hd$p_0<- 1 - pchisq(hd$LRT0, df = 1)
hd$p_2<- 1 - pchisq(hd$LRT2, df = 2)
hd$p_3<- 1 - pchisq(hd$LRT3, df = 1)

hd$bh_0<-p.adjust(hd$p_0,method="BH")
hd$bh_2<-p.adjust(hd$p_2,method="BH")
hd$bh_3<-p.adjust(hd$p_3,method="BH")

hd$lnp0<-log(hd$p_0,base=10)*-1
hd$lnp3<-log(hd$p_3,base=10)*-1

dim(hd[hd$bh_0<=0.05,])#87
dim(hd[hd$bh_2<=0.05,])#0
dim(hd[hd$bh_3<=0.05,])#19

#compare these to those in sca
aj.out1$CompLoc<-paste(aj.out1$Chrom,aj.out1$Pos,sep=".")
fm.out1$CompLoc<-paste(fm.out1$Chrom,fm.out1$Pos,sep=".")
mo.out1$CompLoc<-paste(mo.out1$Chrom,mo.out1$Pos,sep=".")

aj.prune$CompLoc<-paste(aj.prune$Chrom,aj.prune$Pos,sep=".")
fm.prune$CompLoc<-paste(fm.prune$Chrom,fm.prune$Pos,sep=".")
mo.prune$CompLoc<-paste(mo.prune$Chrom,mo.prune$Pos,sep=".")
gw.sum$CompLoc<-paste(gw.sum$Chrom,gw.sum$Pos,sep=".")

fm.both.out<-fm.out1[fm.out1$CompLoc %in% hd[hd$bh_0<=0.05,"Locus"],"CompLoc"]#22
mo.both.out<-mo.out1[mo.out1$CompLoc %in% hd[hd$bh_3<=0.05,"Locus"],"CompLoc"]#5

#############################################################################

####################################PLOT#####################################
#Includes scaffolds
#png("../fst.selection.episodes_redo_all.png",height=300,width=300,
#	units="mm",res=300)
pdf("../fst.selection.episodes_redo_all.pdf",height=11.5,width=11.5)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(1,1,1,0),mgp=c(3,0.5,0), cex=1.5)
mo<-fst.plot(mo.plot, ci.dat=c(mo.top1,0),fst.name="MOM.FEM", 
	chrom.name="Chrom", axis.size=0,bp.name="Pos",
	sig.col=c("purple3","black"))
points(mo$Pos[mo$LocID %in% shared.out$LocID & mo$MOM.FEM >= mo.top1],
	mo$MOM.FEM[mo$LocID %in% shared.out$LocID & mo$MOM.FEM >= mo.top1],
	col="red",pch=8)
points(mo$Pos[mo$CompLoc %in% mo.both.out],
	mo$MOM.FEM[mo$CompLoc %in% mo.both.out],
	col="purple3",pch=5,cex=1)
axis(2,at=seq(0,0.2,0.05),pos=0,las=1,cex.axis=0.75)
legend("top","Females-Inferred Mothers",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(mo[mo$Chrom ==lgs[i],"Pos"]),y=-0.004,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(mo[mo$Chrom ==lgs[i],"Pos"])
}

fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL",
	 chrom.name="Chrom", axis.size=0,bp.name="Pos",
	sig.col=c("green4","black"))
points(fm$Pos[fm$LocID %in% shared.out$LocID & fm$FEM.MAL >= fm.top1],
	fm$FEM.MAL[fm$LocID %in% shared.out$LocID& fm$FEM.MAL >= fm.top1],
	col="red",pch=8)
points(fm$Pos[fm$CompLoc %in% fm.both.out],
	fm$FEM.MAL[fm$CompLoc %in% fm.both.out],
	col="green4",pch=5,cex=1)
axis(2,at=seq(0,0.2,0.1),pos=0,las=1,cex.axis=0.75)
legend("top","Male-Female",bty='n',cex=0.75,text.font=2)

last<-0
for(i in 1:length(lgs)){
	text(x=mean(fm[fm$Chrom ==lgs[i],"Pos"]),y=-0.004,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fm[fm$Chrom ==lgs[i],"Pos"])
}

aj<-fst.plot(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0, bp.name="Pos",
	sig.col=c("dodgerblue","black"))
axis(2,at=c(0,0.025,0.05),pos=0,las=1,cex.axis=0.75)
points(aj$Pos[aj$LocID %in% shared.out$LocID& aj$ADULT.JUVIE >= aj.top1],
	aj$ADULT.JUVIE[aj$LocID %in% shared.out$LocID & 
	aj$ADULT.JUVIE >= aj.top1],
	col="red",pch=8)
legend("top","Adult-Offspring",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(aj[aj$Chrom ==lgs[i],"Pos"]),y=-0.002,
		labels=lgn[i],  adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(aj[aj$Chrom ==lgs[i],"Pos"])
}
mtext("Linkage Group", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","dodgerblue","red"),pch=c(19,19,19,8),
	c("Viability Selection","Sexual Selection","Overall Selection",
		"Shared Outlier"),
	bg="white",ncol=4,box.lty=0)
dev.off()

#png("../fst.selection.episodes_redo_lgs.png",height=300,width=300,units="mm",
	#res=300)
pdf("../fst.selection.episodes_redo_lgs.pdf",height=11.5,width=11.5)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(1,1,1,0),mgp=c(3,0.5,0), cex=1.5)
mo<-fst.plot(mo.plot, ci.dat=c(mo.top1,0),fst.name="MOM.FEM", 
	chrom.name="Chrom", axis.size=0,bp.name="Pos",
	sig.col=c("purple3","black"),groups=lgs)
points(mo$Pos[mo$LocID %in% shared.out$LocID & mo$MOM.FEM >= mo.top1],
	mo$MOM.FEM[mo$LocID %in% shared.out$LocID & mo$MOM.FEM >= mo.top1],
	col="red",pch=8)
points(mo$Pos[mo$CompLoc %in% mo.both.out[substr(mo.both.out,1,2)=="LG"]],
	mo$MOM.FEM[mo$CompLoc %in% mo.both.out[substr(mo.both.out,1,2)=="LG"]],
	col="purple3",pch=5,cex=1)
axis(2,at=seq(0,0.2,0.05),pos=0,las=1,cex.axis=0.75)
legend("top","Females-Inferred Mothers",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(mo[mo$Chrom ==lgs[i],"Pos"]),y=-0.004,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(mo[mo$Chrom ==lgs[i],"Pos"])
}

fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", 
	chrom.name="Chrom", axis.size=0,bp.name="Pos",
	sig.col=c("green4","black"),groups=lgs)
points(fm$Pos[fm$LocID %in% shared.out$LocID & fm$FEM.MAL >= fm.top1],
	fm$FEM.MAL[fm$LocID %in% shared.out$LocID& fm$FEM.MAL >= fm.top1],
	col="red",pch=8)
points(fm$Pos[fm$CompLoc %in% fm.both.out[substr(fm.both.out,1,2)=="LG"]],
	fm$FEM.MAL[fm$CompLoc %in% fm.both.out[substr(fm.both.out,1,2)=="LG"]],
	col="green4",pch=5,cex=1)
axis(2,at=seq(0,0.2,0.1),pos=0,las=1,cex.axis=0.75)
legend("top","Male-Female",bty='n',cex=0.75,text.font=2)

last<-0
for(i in 1:length(lgs)){
	text(x=mean(fm[fm$Chrom ==lgs[i],"Pos"]),y=-0.004,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(fm[fm$Chrom ==lgs[i],"Pos"])
}

aj<-fst.plot(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0, bp.name="Pos",
	sig.col=c("dodgerblue","black"),groups=lgs)
axis(2,at=c(0,0.025,0.05),pos=0,las=1,cex.axis=0.75)
points(aj$Pos[aj$LocID %in% shared.out$LocID& aj$ADULT.JUVIE >= aj.top1],
	aj$ADULT.JUVIE[aj$LocID %in% shared.out$LocID& aj$ADULT.JUVIE >= aj.top1],
	col="red",pch=8)
legend("top","Adult-Offspring",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(aj[aj$Chrom ==lgs[i],"Pos"]),y=-0.002,
		labels=lgn[i],  adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(aj[aj$Chrom ==lgs[i],"Pos"])
}
mtext("Linkage Group", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","dodgerblue","red","black"),
	pch=c(19,19,19,8,5),
	c("Viability Selection","Sexual Selection","Overall Selection",
		"Shared Outlier","Shared with LRT"),
	bg="white",ncol=5,box.lty=0)
dev.off()


png("../male-female.png",height=100,width=300,units="mm",res=300)
par(oma=c(2,2,2,2),mar=c(4,0,0,0))
fm<-fst.plot(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", 
	chrom.name="Chrom", axis.size=0,bp.name="Pos",
	sig.col=c("green4","black"))
axis(2,at=c(0,0.1,0.2,0.3),pos=0,las=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0,line=1)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(fm[fm$Chrom ==lgs[i],"Pos"]),y=-0.002,
		labels=lgs[i], srt=45, adj=1, xpd=TRUE)
	last<-max(fm[fm$Chrom ==lgs[i],"Pos"])
}
dev.off()

##########KELLY/MONNAHAN
pdf("../kelly_analysis_lgs.pdf",height=7.66666,width=11.5)
par(mfrow=c(2,1),oma=c(1,1,0,0),mar=c(1,1,1,0),mgp=c(3,0.5,0), cex=1.5)
hd3<-fst.plot(hd, ci.dat=c(100,-100),fst.name="lnp3", chrom.name="chrom"
	, axis.size=0,bp.name="pos",sig.col=c("green4","black"),groups=lgs)
points(hd3[hd3$bh_3<=0.05 & hd3$chrom %in% lgs,"pos"],
	hd3[hd3$bh_3<=0.05 & hd3$chrom %in% lgs,"lnp3"],
	col="green4",pch=19,cex=0.75)
points(hd3[hd3$bh_3<=0.05 & hd3$Locus %in% mo.out1$CompLoc & 
		hd3$chrom %in% lgs,"pos"],
	hd3[hd3$bh_3<=0.05 & hd3$Locus %in% mo.out1$CompLoc & 
		hd3$chrom %in% lgs,"lnp3"],
	col="green4",pch=5,cex=1)
axis(2,at=seq(0,15,5),pos=0,las=1,cex.axis=0.75)
legend("top","Successful Females-Adults",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(hd3[hd3$chrom ==lgs[i],"pos"]),y=-0.5,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(hd3[hd3$chrom ==lgs[i],"pos"])
}
hd0<-fst.plot(hd, ci.dat=c(100,-100),fst.name="lnp0", chrom.name="chrom"
	, axis.size=0,bp.name="pos",sig.col=c("purple3","black"),groups=lgs)
points(hd0[hd0$bh_0<=0.05 & hd0$chrom %in% lgs,"pos"],
	hd0[hd0$bh_0<=0.05 & hd0$chrom %in% lgs,"lnp0"],
	col="purple3",pch=19,cex=0.75)
points(hd0[hd0$bh_0<=0.05 & hd0$Locus %in% fm.out1$CompLoc & 
		hd0$chrom %in% lgs,"pos"],
	hd0[hd0$bh_0<=0.05 & hd0$Locus %in% fm.out1$CompLoc & 
		hd0$chrom %in% lgs,"lnp0"],
	col="purple3",pch=5,cex=1)
axis(2,at=seq(0,15,5),pos=0,las=1,cex.axis=0.75)
legend("top","Males-Females",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(hd0[hd0$chrom ==lgs[i],"pos"]),y=-0.5,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(hd0[hd0$chrom ==lgs[i],"pos"])
}
mtext("Linkage Group", 1, outer=T, cex=1)
mtext(expression(-log[10]~italic(p)), 2, outer=T, cex=1,las=0)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","black"),pch=c(19,19,5),
	c("Viability Selection","Sexual Selection",
		expression(Shared~with~italic(F)[ST])),
	bg="white",ncol=3,box.lty=0)
dev.off()

pdf("../kelly_analysis_scaffolds.pdf",height=7.66666,width=11.5)
par(mfrow=c(2,1),oma=c(1,1,0,0),mar=c(1,1,1,0),mgp=c(3,0.5,0), cex=1.5)
hd3<-fst.plot(hd, ci.dat=c(100,-100),fst.name="lnp3", chrom.name="chrom"
	, axis.size=0,bp.name="pos",sig.col=c("green4","black"))
points(hd3[hd3$bh_3<=0.05,"pos"],
	hd3[hd3$bh_3<=0.05,"lnp3"],
	col="green4",pch=19,cex=0.75)
points(hd3[hd3$bh_3<=0.05 & hd3$Locus %in% mo.out1$CompLoc,"pos"],
	hd3[hd3$bh_3<=0.05 & hd3$Locus %in% mo.out1$CompLoc,"lnp3"],
	col="green4",pch=5,cex=1)
axis(2,at=seq(0,15,5),pos=0,las=1,cex.axis=0.75)
legend("top","Successful Females-Adults",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(hd3[hd3$chrom ==lgs[i],"pos"]),y=-0.5,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(hd3[hd3$chrom ==lgs[i],"pos"])
}
hd0<-fst.plot(hd, ci.dat=c(100,-100),fst.name="lnp0", chrom.name="chrom"
	, axis.size=0,bp.name="pos",sig.col=c("purple3","black"))
points(hd0[hd0$bh_0<=0.05,"pos"],
	hd0[hd0$bh_0<=0.05,"lnp0"],
	col="purple3",pch=19,cex=0.75)
points(hd0[hd0$bh_0<=0.05 & hd0$Locus %in% fm.out1$CompLoc,"pos"],
	hd0[hd0$bh_0<=0.05 & hd0$Locus %in% fm.out1$CompLoc,"lnp0"],
	col="purple3",pch=5,cex=1)
axis(2,at=seq(0,15,5),pos=0,las=1,cex.axis=0.75)
legend("top","Males-Females",bty='n',cex=0.75,text.font=2)
last<-0
for(i in 1:length(lgs)){
	text(x=mean(hd0[hd0$chrom ==lgs[i],"pos"]),y=-0.5,
		labels=lgn[i], adj=1, xpd=TRUE,srt=90,cex=0.75)
	last<-max(hd0[hd0$chrom ==lgs[i],"pos"])
}
mtext("Linkage Group", 1, outer=T, cex=1)
mtext(expression(-log[10]~italic(p)), 2, outer=T, cex=1,las=0)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","black"),pch=c(19,19,5),
	c("Viability Selection","Sexual Selection",
		expression(Shared~with~italic(F)[ST])),
	bg="white",ncol=3,box.lty=0)
dev.off()


#############################################################################

#################################WRITE TO FILE##############################
write.table(levels(as.factor(aj.unique$LocID)),
	"../biallelic_outliers/AJ_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(fm.unique$LocID)),
	"../biallelic_outliers/FM_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(mo.unique$LocID)),
	"../biallelic_outliers/MO_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(shared.out$LocID)),
	"../biallelic_outliers/shared1.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(aj.mo$LocID)),
	"../biallelic_outliers/ajmo.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(aj.fm$LocID)),
	"../biallelic_outliers/aj.fm.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(fm.mo$LocID)),
	"../biallelic_outliers/fmmo.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
##Monnahan
write.table(fm.out1[fm.out1$CompLoc %in% hd[hd$bh_0<=0.05,"Locus"],
	c("Chrom","Pos","LocID","Locus")],"FM_LRTandFST.txt",sep='\t',
	col.names=T,row.names=F,quote=F)
write.table(mo.out1[mo.out1$CompLoc %in% hd[hd$bh_3<=0.05,"Locus"],
	c("Chrom","Pos","LocID","Locus")],"MO_LRTandFST.txt",sep='\t',
	col.names=T,row.names=F,quote=F)

###WRITE SCAFFOLD, START, and END to file 
aj.unique<-aj.plot[aj.plot$Locus %in% aj.unique$Locus,]
aj.rad.region<-data.frame(aj.unique$Chrom,as.numeric(aj.unique$Pos)-2500,
	as.numeric(aj.unique$Pos)+2500)
write.table(aj.rad.region,"../biallelic_outliers/rad_region/aj_extract.txt",
	quote=F,col.names=F,row.names=F,sep="\t")
fm.unique<-fm.plot[fm.plot$Locus %in% fm.unique$Locus,]
fm.rad.region<-data.frame(fm.unique$Chrom,as.numeric(fm.unique$Pos)-2500,
	as.numeric(fm.unique$Pos)+2500)
write.table(fm.rad.region,"../biallelic_outliers/rad_region/fm_1extract_redo.txt",
	quote=F,col.names=F,row.names=F,sep="\t")
mo.unique<-mo.plot[mo.plot$Locus %in% mo.unique$Locus,]
mo.rad.region<-data.frame(mo.unique$Chrom,as.numeric(mo.unique$Pos)-2500,
	as.numeric(mo.unique$Pos)+2500)
write.table(mo.rad.region,"../biallelic_outliers/rad_region/mo_1extract.txt",
	quote=F,col.names=F,row.names=F,sep="\t")

aj.mo<-aj.plot[aj.plot$Locus %in% aj.mo$Locus,]
aj.mo.region<-data.frame(aj.mo$Chrom, as.numeric(aj.mo$Pos)-2500,
	as.numeric(aj.mo$Pos)+2500)
write.table(aj.mo.region,"../biallelic_outliers/rad_region/ajmo_extract.txt",
	quote=F,col.names=F,row.names=F,sep='\t')
fm.mo<-fm.plot[fm.plot$Locus %in% fm.mo$Locus,]
fm.mo.region<-data.frame(fm.mo$Chrom, as.numeric(fm.mo$Pos)-2500,
	as.numeric(fm.mo$Pos)+2500)
write.table(fm.mo.region,"../biallelic_outliers/rad_region/fmmo_extract.txt",
	quote=F,col.names=F,row.names=F,sep='\t')
aj.fm<-aj.plot[aj.plot$Locus %in% aj.fm$Locus,]
aj.fm.region<-data.frame(aj.fm$Chrom, as.numeric(aj.fm$Pos)-2500,
	as.numeric(aj.fm$Pos)+2500)
write.table(aj.fm.region,"../biallelic_outliers/rad_region/ajfm_extract.txt",
	quote=F,col.names=F,row.names=F,sep='\t')

sharedregion<-aj.plot[aj.plot$LocID %in% shared.out$LocID,]
shared.region<-data.frame(sharedregion$Chrom,as.numeric(sharedregion$Pos-2500),
	as.numeric(sharedregion$Pos)+2500)
write.table(shared.region,"../biallelic_outliers/rad_region/shared_2extract.txt",
	quote=F,col.names=F,row.names=F,sep='\t')


##Monnahan
lrt.fm.region<-data.frame(hd[hd$bh_0<=0.05,"chrom"],
	as.numeric(hd[hd$bh_0<=0.05,"pos"])-2500,
	as.numeric(hd[hd$bh_0<=0.05,"pos"])+2500)
write.table(lrt.fm.region,"../biallelic_outliers/rad_region/lrt_fm.txt",
	quote=F,col.names=F,row.names=F,eol='\n')
lrt.mo.region<-data.frame(hd[hd$bh_3<=0.05,"chrom"],
	as.numeric(hd[hd$bh_3<=0.05,"pos"])-2500,
	as.numeric(hd[hd$bh_3<=0.05,"pos"])+2500)
write.table(lrt.mo.region,"../biallelic_outliers/rad_region/lrt_mo.txt",
	quote=F,col.names=F,row.names=F,eol='\n')
#############################################################################

########################LOOK INTO THE EXTREME OUTLIERS#######################
fm.extreme<-fm[fm$FEM.MAL >=0.11,c("Locus","Chrom","Pos","LocID","FEM.MAL")]
mo.extreme<-mo[mo$MOM.FEM >= 0.08,c("Locus","Chrom","Pos","LocID","MOM.FEM")]
aj.extreme<-aj[aj$ADULT.JUVIE >= 0.02,c("Locus","Chrom","Pos","LocID",
	"ADULT.JUVIE")]

mo.ex.sum<-gw.sum[gw.sum$Locus %in% mo.extreme$Locus & 
	gw.sum$Pop %in% c("FEM","MOM"),]

fm.ex.sum<-gw.sum[gw.sum$Locus %in% fm.extreme$Locus & 
	gw.sum$Pop %in% c("FEM","MAL"),]

aj.ex.sum<-gw.sum[gw.sum$Locus %in% aj.extreme$Locus &
	gw.sum$Pop %in% c("ADULT","JUVIE"),]

par(mfrow=c(3,2))
hist(aj.ex.sum$Allele1Freq[aj.ex.sum$Pop=="ADULT"],xlab="",	ylab="",
	main="Adults")
hist(aj.ex.sum$Allele1Freq[aj.ex.sum$Pop=="JUVIE"],xlab="",	ylab="",
	main="Offspring")
hist(mo.ex.sum$Allele1Freq[mo.ex.sum$Pop=="MOM"],xlab="",
	ylab="",main="Mothers")
hist(mo.ex.sum$Allele1Freq[mo.ex.sum$Pop=="FEM"],xlab="",
	ylab="",main="Females")
hist(fm.ex.sum$Allele1Freq[fm.ex.sum$Pop=="JUVIE"],xlab="",	ylab="",
	main="Offspring")
hist(fm.ex.sum$Allele1Freq[fm.ex.sum$Pop=="MAL"],xlab="",
	ylab="",main="Males")

col.order<-c("LocusID","Chrom","BP","NumSNPs","NumSNPsOut","Pop1N","Pop2N",
	"Pop1AF","Pop2AF","Pop1Ho","Pop2Ho","Pop1Hs","Pop2Hs","SNP1","SNP2",
	"SNP3","SNP4","Seq")

aj.ex.dat<-extract.info(aj.ex.sum,tags,snps)
aj.ex.dat<-aj.ex.dat[,col.order]
aj.ex.dat$Comparison<-"Adult-Offspring"
mo.ex.dat<-extract.info(mo.ex.sum,tags,snps)
mo.ex.dat<-mo.ex.dat[,col.order]
mo.ex.dat$Comparison<-"Mothers-Females"
fm.ex.dat<-extract.info(fm.ex.sum,tags,snps)
fm.ex.dat<-fm.ex.dat[,col.order]
fm.ex.dat$Comparison<-"Female-Male"

extreme.outliers<-rbind(aj.ex.dat,mo.ex.dat,fm.ex.dat)
write.csv(extreme.outliers,"ExtremeOutliers.csv",row.names=F)

#####All Outliers
aj.null<-aj.prune[!(aj.prune$Locus %in% aj.out$Locus),]
aj.null.sum<-gw.sum[gw.sum$Locus %in% aj.null$Locus,]
aj.null.sum<-aj.null.sum[aj.null.sum$Pop %in% c("ADULT","JUVIE"),]
aj.all<-extract.info(aj.null.sum,tags,snps)

fm.null<-fm.prune[!(fm.prune$Locus %in% fm.out$Locus),]
fm.null.sum<-gw.sum[gw.sum$Locus %in% fm.null$Locus,]
fm.null.sum<-fm.null.sum[fm.null.sum$Pop %in% c("MAL","FEM"),]
fm.all<-extract.info(fm.null.sum,tags,snps)

mo.null<-mo.prune[!(mo.prune$Locus %in% mo.out$Locus),]
mo.null.sum<-gw.sum[gw.sum$Locus %in% mo.null$Locus,]
mo.null.sum<-mo.null.sum[mo.null.sum$Pop %in% c("MOM","FEM"),]
mo.all<-extract.info(mo.null.sum,tags,snps)

aj.out.dat<-gw.sum[gw.sum$Locus %in% aj.out$Locus,]
aj.out.dat<-aj.out.dat[aj.out.dat$Pop %in% c("ADULT","JUVIE"),]
aj.out.dat<-extract.info(aj.out.dat,tags,snps)
mo.out.dat<-gw.sum[gw.sum$Locus %in% mo.out$Locus,]
mo.out.dat<-mo.out.dat[mo.out.dat$Pop %in% c("MOM","FEM"),]
mo.out.dat<-extract.info(mo.out.dat,tags,snps)
fm.out.dat<-gw.sum[gw.sum$Locus %in% fm.out$Locus,]
fm.out.dat<-fm.out.dat[fm.out.dat$Pop %in% c("MAL","FEM"),]
fm.out.dat<-extract.info(fm.out.dat,tags,snps)
nrad<-c(nrow(aj.null),nrow(aj.out.dat),nrow(aj.ex.dat),
	nrow(mo.all),nrow(mo.out.dat),nrow(mo.ex.dat),
	nrow(fm.all),nrow(fm.out.dat),nrow(fm.ex.dat))
library(scales)
png("SNPsPerRADLocus.png", height=7,width=7,units="in",res=300)
boxplot(aj.all$NumSNPs,aj.out.dat$NumSNPs,aj.ex.dat$NumSNPs,
	mo.all$NumSNPs,mo.out.dat$NumSNPs,mo.ex.dat$NumSNPs,
	fm.all$NumSNPs,fm.out.dat$NumSNPs,fm.ex.dat$NumSNPs,
	col=c("grey",alpha("green4",0.5),"green4",
		"grey",alpha("purple3",0.5),"purple3",
		"grey",alpha("dodgerblue",0.5),"dodgerblue"),
	names=rep("",9),ylim=c(0,70),axes=F)
axis(2,at=seq(0,70,10),las=1)
mtext("Number of SNPs Per RAD tag",2,line=2)
text(1:9,par("usr")[3] - 1, srt = 45, adj = 1,xpd = TRUE,
     labels = c("Adult-Off Null","Adult-Off Outliers","Adult-Off Extreme", 
	"Fem-Mom Null","Fem-Mom Outliers","Fem-Mom Extreme",
	"Mal-Fem Null","Mal-Fem Outliers","Mal-Fem Extreme"))
text(1:9,y=66,labels=nrad)
dev.off()

snps.per.rad<-data.frame(rbind(cbind(aj.all$NumSNPs,"Adult-Off","Null"),
	cbind(aj.out.dat$NumSNPs,"Adult-Off","Outliers"),
	cbind(aj.ex.dat$NumSNPs,"Adul-Off","Extreme"),
	cbind(mo.all$NumSNPs,"Fem-Mom","Null"),
	cbind(mo.out.dat$NumSNPs,"Fem-Mom","Outliers"),
	cbind(mo.ex.dat$NumSNPs,"Fem-Mom","Extreme"),
	cbind(fm.all$NumSNPs,"Mal-Fem","Null"),
	cbind(fm.out.dat$NumSNPs,"Mal-Fem","Outliers"),
	cbind(fm.ex.dat$NumSNPs,"Mal-Fem","Extreme")))
colnames(snps.per.rad)<-c("NumSNPs","Comparison","SNPType")
spr.lm<-lm(as.numeric(NumSNPs)~Comparison+SNPType,data=snps.per.rad)
#> anova(spr.lm)
#Analysis of Variance Table
#
#Response: as.numeric(NumSNPs)
#              Df   Sum Sq Mean Sq F value  Pr(>F)    
#Comparison     3     1891   630.4   3.211 0.02197 *  
#SNPType        2    15494  7746.9  39.456 < 2e-16 ***
#Residuals  57949 11377719   196.3   

#Do they have more Ns?
Ns<-data.frame(rbind(do.call(cbind,lapply(
	apply(aj.all[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(aj.out.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(aj.ex.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(mo.all[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(mo.out.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(mo.ex.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(fm.all[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(fm.out.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(fm.ex.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length))))
Ns$Comparison<-c(rep("Adult-Off",3),rep("Fem-Mom",3),rep("Mal-Fem",3))
Ns$SNPType<-rep(c("Null","Outliers","Extreme"),3)
Ns$NumLoci<-c(nrow(aj.all),nrow(aj.out.dat),nrow(aj.ex.dat),
	nrow(mo.all),nrow(mo.out.dat),nrow(mo.ex.dat),
	nrow(fm.all),nrow(fm.out.dat),nrow(fm.ex.dat))
Ns$PropN<-rowSums(Ns[,1:4])/Ns$NumLoci

png("PropN.png",height=7,width=7,units="in",res=300)
bp<-barplot(Ns$PropN,
	col=c("grey",alpha("green4",0.5),"green4",
		"grey",alpha("purple3",0.5),"purple3",
		"grey",alpha("dodgerblue",0.5),"dodgerblue"),
	names="",ylim=c(0,1),las=1,ylab="Number of SNPs Per RAD tag")
text(bp,par("usr")[3] - 0.01, srt = 45, adj = 1,xpd = TRUE,
     labels = c("Adult-Off Null","Adult-Off Outliers","Adult-Off Extreme", 
		"Fem-Mom Null","Fem-Mom Outliers","Fem-Mom Extreme",
		"Mal-Fem Null","Mal-Fem Outliers","Mal-Fem Extreme"))
text(bp,y=0.05,labels=Ns$NumLoci,srt=90)
dev.off()

#For outliers on loci with a bunch of SNPs, are most of the SNPs significant?
hist(aj.out.dat$NumSNPsOut/aj.out.dat$NumSNPs)
hist(fm.out.dat$NumSNPsOut/fm.out.dat$NumSNPs)
hist(mo.out.dat$NumSNPsOut/mo.out.dat$NumSNPs)

#OUTLIERS
hd.fm<-hd[hd$bh_0<=0.05,]#87
hd.mo<-hd[hd$bh_3<=0.05,]#19

##MATCH TO OTHERS
vcf$CompLoc<-paste(vcf$`#CHROM`,vcf$POS,sep=".")
vcf$Locus<-paste(vcf$`#CHROM`,vcf$ID,vcf$POS,sep=".")
info<-vcf[,c("#CHROM","POS","ID","CompLoc","Locus")]

snp.dat<-parse.snps(snps)
loc.dat<-merge(info,snp.dat,by.x="ID",by.y="LocusID",all=T)
loc.dat<-merge(loc.dat,tags,by.x="ID",by.y="LocusID",all=T)
#MATCH TO LOCUS INFO
hd.locus<-hd[hd$Locus %in% loc.dat$CompLoc,]
hd.fm.dat<-loc.dat[loc.dat$CompLoc %in% hd.fm$Locus,]#83
hd.mo.dat<-loc.dat[loc.dat$CompLoc %in% hd.mo$Locus,]#19

fm.lrt.dat<-data.frame(Chrom=hd.fm.dat$Chrom,ID=hd.fm.dat$ID,
	Pos=hd.fm.dat$POS,LocusID=hd.fm.dat$Locus,
	NumSNPs=hd.fm.dat$NumSNPs,Seq=hd.fm.dat$Seq)
write.table(fm.lrt.dat,"../biallelic_outliers/fm_lrt_dat.txt",col.names=T,
	quote=F)
mo.lrt.dat<-data.frame(Chrom=hd.mo.dat$Chrom,ID=hd.mo.dat$ID,
	Pos=hd.mo.dat$POS,LocusID=hd.mo.dat$Locus,
	NumSNPs=hd.mo.dat$NumSNPs,Seq=hd.mo.dat$Seq)
write.table(mo.lrt.dat,"../biallelic_outliers/mo_lrt_dat.txt",col.names=T,
	quote=F)
#Next up: Merge with Blast2Go results.

#############################################################################

#############################BLAST2GO ANNOTATIONS############################
setwd("../biallelic_outliers/blast")
out.cnames<-c("LocID","NumSNPs","Chrom.x","Pos","Description","Length",
	"X.Hits","e.Value","sim.mean","X.GO","GO.Names.list","Enzyme.Codes.list","Seq")
out.names<-c("LocID","NumSNPs","Chrom","Pos","Description","Length",
	"NumHits","e.Value","sim.mean","NumGO","GO.Names.list","Enzyme.Codes.list","Seq")
blast2go.files<-list.files(pattern="blast2go")
ajfm.b2g<-read.delim("ajfm_blast2go.txt",sep='\t',header=T)
ajfm.b2g$Comparison<-"AO-FM"
ajmo.b2g<-read.delim("ajmo_blast2go.txt",sep='\t',header=T)
ajmo.b2g$Comparison<-"AO-MO"
fmmo.b2g<-read.delim("fmmo_blast2go.txt",sep='\t',header=T)
fmmo.b2g$Comparison<-"MO-FM"
shared.b2g<-read.delim("sharedall_blast2go.txt",sep='\t',header=T,stringsAsFactors=F)
shared1.b2g<-read.delim("shared_blast2go.txt",sep='\t',header=T,stringsAsFactors=F)
shared.b2g<-rbind(shared.b2g,shared1.b2g)


shared.out$SeqName<-"LG21_1363556-1368556"
shared.out$Comparison<-"ALL"

aj.fm$start<-aj.fm$Pos-2500
aj.fm$end<-aj.fm$Pos+2500
aj.fm$start[aj.fm$start<0]<-0
aj.fm$SeqName<-paste(aj.fm$Chrom,"_",aj.fm$start,"-",aj.fm$end,sep="")
#ajfm<-merge(aj.fm,ajfm.b2g,by="SeqName")
#ajfm<-merge(ajfm,aj.out.dat,by.x="LocID",by.y="LocusID")
#ajfm<-ajfm[,out.cnames]
#colnames(ajfm)<-out.names
aj.fm$Comparison<-"AdultOff-MalFem"
aj.fm<-aj.fm[!duplicated(aj.fm$LocID),]

aj.mo$start<-aj.mo$Pos-2500
aj.mo$end<-aj.mo$Pos+2500
aj.mo$start[aj.mo$start<0]<-0
aj.mo$SeqName<-paste(aj.mo$Chrom,"_",aj.mo$start,"-",aj.mo$end,sep="")
#ajmo<-merge(aj.mo,ajmo.b2g,by="SeqName")
#ajmo<-merge(ajmo,mo.out.dat,by.x="LocID",by.y="LocusID")
#ajmo<-ajmo[,out.cnames]
#colnames(ajmo)<-out.names
aj.mo$Comparison<-"AdultOff-FemMom"
aj.mo<-aj.mo[!duplicated(aj.mo$LocID),]

fm.mo$start<-fm.mo$Pos-2500
fm.mo$end<-fm.mo$Pos+2500
fm.mo$start[fm.mo$start<0]<-0
fm.mo$SeqName<-paste(fm.mo$Chrom,"_",fm.mo$start,"-",fm.mo$end,sep="")
fm.mo$Comparison<-"FemMal-FemMom"
fm.mo<-fm.mo[!duplicated(fm.mo$LocID),]

shared<-rbind(aj.fm,aj.mo,fm.mo,shared.out)
shared.blast<-merge(shared,shared.b2g,by="SeqName")
all.out.dat<-rbind(mo.out.dat,aj.out.dat,fm.out.dat)
all.out.dat<-all.out.dat[!duplicated(all.out.dat$LocusID),]
shared.blast<-merge(shared.blast,all.out.dat,by.x="LocID",by.y="LocusID")
shared.blast<-cbind(shared.blast[,"Comparison"],shared.blast[,out.cnames])
colnames(shared.blast)<-c("Comparison",out.names)

write.table(shared.blast,"../../S1_SharedOutliersBlast.txt",sep='\t',col.names=T,
	row.names=F,quote=F)

aj.unique<-aj.plot[aj.plot$Locus %in% aj.unique$Locus,
	c("Chrom","Pos","LocID","Locus")]
aj.unique$start<-aj.unique$Pos-2500
aj.unique$end<-aj.unique$Pos+2500
aj.unique$start[aj.unique$start < 0]<-"0"
aj.unique$SeqName<-paste(aj.unique$Chrom,"_",aj.unique$start,"-",
	aj.unique$end,sep="")
aj.unique<-aj.unique[!duplicated(aj.unique$LocID),]
aj.b2g<-read.delim("aj_blast2go.txt",header=T,sep='\t')
aj.blast<-merge(aj.unique,aj.b2g,by="SeqName")
aj.blast<-merge(aj.blast,aj.out.dat,by.x="LocID",by.y="LocusID")
aj.blast<-aj.blast[,out.cnames]
colnames(aj.blast)<-out.names
aj.blast$Comparison<-"Adults-Offspring"

fm.unique<-fm.plot[fm.plot$Locus %in% fm.unique$Locus,
	c("Chrom","Pos","LocID","Locus")]
fm.unique$start<-fm.unique$Pos-2500
fm.unique$end<-fm.unique$Pos+2500
fm.unique$start[fm.unique$start < 0]<-"0"
fm.unique$SeqName<-paste(fm.unique$Chrom,"_",fm.unique$start,"-",
	fm.unique$end,sep="")
fm.unique<-fm.unique[!duplicated(fm.unique$LocID),]
fm.b2g<-read.delim("fm_blast2go.txt",header=T,sep='\t')
fm.blast<-merge(fm.unique,fm.b2g,by="SeqName")
fm.blast<-merge(fm.blast,fm.out.dat,by.x="LocID",by.y="LocusID")
fm.blast<-fm.blast[,out.cnames]
colnames(fm.blast)<-out.names
fm.blast$Comparison<-"Males-Females"

mo.unique<-mo.plot[mo.plot$Locus %in% mo.unique$Locus,
	c("Chrom","Pos","LocID","Locus")]
mo.unique$start<-mo.unique$Pos-2500
mo.unique$end<-mo.unique$Pos+2500
mo.unique$start[mo.unique$start < 0]<-"0"
mo.unique$SeqName<-paste(mo.unique$Chrom,"_",mo.unique$start,"-",
	mo.unique$end,sep="")
mo.unique<-mo.unique[!duplicated(mo.unique$LocID),]
mob2g<-read.delim("mo_blast2go.txt",header=T,sep='\t')
mo.b2g<-merge(mo.unique,mob2g,by="SeqName")
mo.blast<-merge(mo.b2g,mo.out.dat,by.x="LocID",by.y="LocusID")
mo.blast<-mo.blast[,out.cnames]
colnames(mo.blast)<-out.names
mo.blast$Comparison<-"Mothers-Females"

unique<-rbind(aj.blast,fm.blast,mo.blast)
write.table(unique,"../../S2_UniqueOutliersBlast.txt",row.names=F,col.names=T,
	sep='\t',quote=F)
#############################################################################

#######################COMPARE TO PSTFST SIGNIFICANT LOCI####################
pstfst<-read.csv(row.names=1,
	"E:/ubuntushare/popgen/sw_results/pstfst/pstfst_loci_summary.csv")
pstfst.bands<-pstfst[pstfst$Trait=="Bands",]
length(pstfst.bands$scaff[pstfst.bands$scaff %in% aj.out$Chrom])
length(pstfst.bands$scaff[pstfst.bands$scaff %in% mo.out$Chrom])
length(pstfst.bands$scaff[pstfst.bands$scaff %in% bj.out$Chrom])
pstfst.svl<-pstfst[pstfst$Trait=="SVL",]
length(pstfst.svl$scaff[pstfst.svl$scaff %in% aj.out$Chrom])
length(pstfst.svl$scaff[pstfst.svl$scaff %in% mo.out$Chrom])
length(pstfst.svl$scaff[pstfst.svl$scaff %in% bj.out$Chrom])
popgen.out<-read.table(header=T,
	"E:/ubuntushare/popgen/sw_results/AllOutliers.txt")
length(popgen.out$scaffold[popgen.out$scaffold %in% aj.out$Chrom])
length(popgen.out$scaffold[popgen.out$scaffold %in% mo.out$Chrom])
length(popgen.out$scaffold[popgen.out$scaffold %in% bj.out$Chrom])
############################################################################


