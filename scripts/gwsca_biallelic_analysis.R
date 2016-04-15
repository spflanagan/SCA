#Author: Sarah P. Flanagan
#Date: 11 February 2016
#Purpose: Analyze the output from gwsca_biallelic.

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")
gw.fst<-read.delim("gwsca_fsts_new.txt")
gw.sum<-read.delim("gwsca_summary_new.txt")
#gw.fst<-gw.fst[,c("Chrom","Pos","ADULT.JUVIE","MAL.FEM","PREGGER.OFF",
#	"MOM.FEM")]
gw.fst$Locus<-paste(gw.fst$Chrom,gw.fst$LocID,gw.fst$Pos,sep=".")
gw.sum$Locus<-paste(gw.sum$Chrom,gw.sum$LocID,gw.sum$Pos,sep=".")

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
sum.prune<-gw.hwe[gw.hwe$Pop=="ADULT" | gw.hwe$Pop=="JUVIE" | 
	gw.hwe$Pop == "POP",]
sum.sum<-tapply(sum.prune$N,sum.prune$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 393]
sum.prune<-gw.hwe[gw.hwe$Locus %in% names(sum.sum),]

##prune based on representation in the groups
sum.list<-split(sum.prune, sum.prune$Pop)
adt.n<-sum.list$ADULT[sum.list$ADULT$N > 216 & !is.na(sum.list$ADULT$Hs),]
juv.n<-sum.list$JUVIE[sum.list$JUVIE$N > 157 & !is.na(sum.list$JUVIE$Hs),]
fem.n<-sum.list$FEM[sum.list$FEM$N > 57& !is.na(sum.list$FEM$Hs),]
mal.n<-sum.list$MAL[sum.list$MAL$N>159& !is.na(sum.list$MAL$Hs),]
mom.n<-sum.list$MOM[sum.list$MOM$N>133& !is.na(sum.list$MOM$Hs),]
brd.n<-sum.list$BREEDER[sum.list$BREEDER$N>290& !is.na(sum.list$BREEDER$Hs),]
#pop.n<-sum.list$POP[sum.list$POP$N>57& !is.na(sum.list$POP$Hs),]
prg.n<-sum.list$PREGGER[sum.list$PREGGER$N>159& !is.na(sum.list$PREGGER$Hs),]
non.n<-sum.list$NONPREG[sum.list$NONPREG$N>14& !is.na(sum.list$NONPREG$Hs),]

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]
brd.n<-brd.n[brd.n$Allele1Freq > 0.05 & brd.n$Allele1Freq < 0.95,]
#pop.n<-pop.n[pop.n$Allele1Freq > 0.05 & pop.n$Allele1Freq < 0.95,]
prg.n<-prg.n[prg.n$Allele1Freq > 0.05 & prg.n$Allele1Freq < 0.95,]
non.n<-non.n[non.n$Allele1Freq > 0.05 & non.n$Allele1Freq < 0.95,]

#write these out to read in for error analysis
write.table(prg.n[,2:4],"LociInPregMales.txt",col.names=T,row.names=F,quote=F)
write.table(juv.n[,2:4],"LociInOffspring.txt",col.names=T,row.names=F,quote=F)
#comparisons
#viability
aj.prune<-gw.fst[gw.fst$Locus %in% adt.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
aj.prune<-aj.prune[aj.prune$ADULT.JUVIE>0,]
#fm.prune<-gw.fst[gw.fst$Locus %in% mal.n$Locus&gw.fst$Locus %in% fem.n$Locus, ]
#fm.prune<-fm.prune[fm.prune$FEM.MAL>0,]
#sexual
mo.prune<-gw.fst[gw.fst$Locus %in% fem.n$Locus&gw.fst$Locus %in% mom.n$Locus, ]
mo.prune<-mo.prune[mo.prune$FEM.MOM>0,]
#gametic
bj.prune<-gw.fst[gw.fst$Locus %in% brd.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
bj.prune<-bj.prune[bj.prune$JUVIE.BREEDER>0,]


write.table(pj.prune, "pj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(fm.prune, "fm.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(mo.prune, "mo.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(aj.prune,"aj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(bj.prune,"bj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')

aj.plot<-aj.prune[order(aj.prune$ADULT.JUVIE),] #ascending
aj.top1<-aj.plot[round(nrow(aj.plot)*0.99),"ADULT.JUVIE"]
aj.out1<-aj.plot[aj.plot$ADULT.JUVIE >= aj.top1,]
#fm.plot<-fm.prune[order(fm.prune$FEM.MAL),]#ascending
#fm.top1<-fm.plot[round(nrow(fm.plot)*0.99),"FEM.MAL"]
#fm.out1<-fm.prune[fm.prune$FEM.MAL >= fm.top1,]
mo.plot<-mo.prune[order(mo.prune$FEM.MOM),]#ascending
mo.top1<-mo.plot[round(nrow(mo.plot)*0.99),"FEM.MOM"]
mo.out1<-mo.prune[mo.prune$FEM.MOM >= mo.top1,]
#pj.plot<-pj.prune[order(pj.prune$JUVIE.PREGGER),] #ascending
#pj.top1<-pj.plot[round(nrow(pj.plot)*0.99),"JUVIE.PREGGER"]
#pj.out1<-pj.prune[pj.prune$JUVIE.PREGGER >= pj.top1,]
bj.plot<-bj.prune[order(bj.prune$JUVIE.BREEDER),]#ascending
bj.top1<-bj.plot[round(nrow(bj.plot)*0.99),"JUVIE.BREEDER"]
bj.out1<-bj.prune[bj.prune$JUVIE.BREEDER >= bj.top1,]

#plot with the top1%
png("fst.top1.comp3_redo.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)

#fm<-plot.fsts(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
#	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
#legend("top","Male-Female", cex=0.75,bty="n")
aj<-plot.fsts(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
legend("top","Adult-Offspring", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)

mo<-plot.fsts(mo.plot, ci.dat=c(mo.top1,0),fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)

#pj<-plot.fsts(pj.plot, ci.dat=c(pj.top1,0),fst.name="JUVIE.PREGGER", chrom.name="Chrom"
#	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
#legend("top","Pregnant Males-Offspring", cex=0.75,bty="n")
bj<-plot.fsts(bj.plot, ci.dat=c(bj.top1,0),fst.name="JUVIE.BREEDER", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
#legend("top","Breeding Adults-Offspring", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()



###COMPARISONS
aj.out<-aj[aj$ADULT.JUVIE >= aj.top1[1],]
#fm.out<-fm[fm$FEM.MAL >= fm.top1[1],]
mo.out<-mo[mo$FEM.MOM >= mo.top1[1],]
#pj.out<-pj[pj$JUVIE.PREGGER >= pj.top1[1],]
bj.out<-bj[bj$JUVIE.BREEDER >= bj.top1[1],]

aj.unique<-aj.out[!(aj.out$Locus %in% bj.out$Locus) & 
	!(aj.out$Locus %in% mo.out$Locus),]
#fm.unique<-fm.out[#!(fm.out$Locus %in% aj.out$Locus) & 
#	!(fm.out$Locus %in% mo.out$Locus) & !(fm.out$Locus %in% pj.out$Locus),]
mo.unique<-mo.out[!(mo.out$Locus %in% aj.out$Locus) &
	!(mo.out$Locus %in% aj.out$Locus),]
#pj.unique<-pj.out[#!(pj.out$Locus %in% aj.out$Locus) &
#	!(pj.out$Locus %in% fm.out$Locus) & !(pj.out$Locus %in% mo.out$Locus),]
bj.unique<-bj.out[!(bj.out$Locus %in% aj.out$Locus) &
	!(bj.out$Locus %in% mo.out$Locus),]
shared<-aj.out[(aj.out$LocID %in% mo.out$LocID) & 
	(aj.out$LocID %in% bj.out$LocID),]

aj.bj<-aj.out[(aj.out$LocID %in% bj.out$LocID),]
aj.mo<-aj.out[(aj.out$LocID %in% mo.out$LocID),]
bj.mo<-bj.out[bj.out$LocID %in% mo.out$LocID,]

png("fst.selection.episodes_redo.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
aj<-plot.fsts(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos",sig.col=c("green4","black"))
points(aj$Pos[aj$LocID %in% shared$LocID& aj$ADULT.JUVIE >= aj.top1],
	aj$ADULT.JUVIE[aj$LocID %in% shared$LocID& aj$ADULT.JUVIE >= aj.top1],
	col="red",pch=8)
legend("top","Adult-Offspring",bty='n',cex=0.75)
#fm<-plot.fsts(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("green4","black"))
#points(vi.out$Pos,vi.out$FEM.MAL,col="dodgerblue",pch=16)
#points(fm.out$Pos[fm.out$Locus %in% shared$Locus],
#	fm.out$FEM.MAL[fm.out$Locus %in% shared$Locus],col="red",pch=8)
#legend("top","Male-Female",bty='n',cex=0.75)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)

mo<-plot.fsts(mo.plot, ci.dat=c(mo.top1,0),fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("purple3","black"))
points(mo$Pos[mo$LocID %in% shared$LocID & mo$FEM.MOM >= mo.top1],
	mo$FEM.MOM[mo$LocID %in% shared$LocID& mo$FEM.MOM >= mo.top1],
	col="red",pch=8)
legend("top","Females-Inferred Mothers",bty='n',cex=0.75)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)

#pj<-plot.fsts(pj.plot, ci.dat=c(pj.top1,0),fst.name="JUVIE.PREGGER", chrom.name="Chrom"
#	, axis.size=0.75,bp.name="Pos",sig.col=c("dodgerblue","black"))
#points(pj.out$Pos[pj.out$Locus %in% shared$Locus],
#	pj.out$JUVIE.PREGGER[pj.out$Locus %in% shared$Locus],col="red",pch=8)
#legend("top","Pregnant Males-Offspring",bty='n',cex=0.75)
bj<-plot.fsts(bj.plot, ci.dat=c(bj.top1,0),fst.name="JUVIE.BREEDER", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("dodgerblue","black"))
points(bj$Pos[bj$LocID %in% shared$LocID & bj$JUVIE.BREEDER >= bj.top1],
	bj$JUVIE.BREEDER[bj$LocID %in% shared$LocID& bj$JUVIE.BREEDER >= bj.top1]
	,col="red",pch=8)
legend("top","Breeding Adults-Offspring",bty='n',cex=0.75)

mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","dodgerblue","red"),pch=c(19,19,19,8),
	c("Viability Selection","Sexual Selection","Gametic Selection",
		"Shared in all"),
	bg="white",ncol=4,box.lty=0)
dev.off()

##WRITE TO FILE
write.table(levels(as.factor(aj.unique$LocID)),
	"../biallelic_outliers/AJ_outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
#write.table(levels(as.factor(fm.unique$LocID)),
#	"../biallelic_outliers/FM_1outliers.txt",quote=F,
#	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(mo.unique$LocID)),
	"../biallelic_outliers/MO_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
#write.table(levels(as.factor(pj.unique$LocID)),
#	"../biallelic_outliers/PJ_1outliers.txt",quote=F,
#	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(bj.out$LocID)),
	"../biallelic_outliers/BJ_outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(shared$LocID)),
	"../biallelic_outliers/shared1.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(aj.mo$LocID)),
	"../biallelic_outliers/ajmo.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(aj.bj$LocID)),
	"../biallelic_outliers/aj.bj.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(bj.mo$LocID)),
	"../biallelic_outliers/bjmo.shared.txt",quote=F,
	col.names=F,row.names=F,sep='\t')

aj.unique<-aj.plot[aj.plot$Locus %in% aj.unique$Locus,]
aj.rad.region<-data.frame(aj.unique$Chrom,as.numeric(aj.unique$Pos)-2500,
	as.numeric(aj.unique$Pos)+2500)
write.table(aj.rad.region,"../biallelic_outliers/rad_region/aj_extract.sh",
	quote=F,col.names=F,row.names=F,sep="\t")
#fm.unique<-fm.plot[fm.plot$Locus %in% fm.unique$Locus,]
#fm.rad.region<-data.frame(fm.unique$Chrom,as.numeric(fm.unique$Pos)-2500,
#	as.numeric(fm.unique$Pos)+2500)
#write.table(fm.rad.region,"../biallelic_outliers/rad_region/fm_1extract_redo.sh",
#	quote=F,col.names=F,row.names=F,sep="\t")
mo.unique<-mo.plot[mo.plot$Locus %in% mo.unique$Locus,]
mo.rad.region<-data.frame(mo.unique$Chrom,as.numeric(mo.unique$Pos)-2500,
	as.numeric(mo.unique$Pos)+2500)
write.table(mo.rad.region,"../biallelic_outliers/rad_region/mo_1extract.sh",
	quote=F,col.names=F,row.names=F,sep="\t")
bj.unique<-bj.plot[bj.plot$Locus %in% bj.unique$Locus,]
bj.rad.region<-data.frame(bj.unique$Chrom,as.numeric(bj.unique$Pos)-2500,
	as.numeric(bj.unique$Pos)+2500)
write.table(bj.rad.region,"../biallelic_outliers/rad_region/bj_1extract.sh",
	quote=F,col.names=F,row.names=F,sep="\t")

aj.mo<-aj.plot[aj.plot$Locus %in% aj.mo$Locus,]
aj.mo.region<-data.frame(aj.mo$Chrom, as.numeric(aj.mo$Pos)-2500,
	as.numeric(aj.mo$Pos)+2500)
write.table(aj.mo.region,"../biallelic_outliers/rad_region/ajmo_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')
bj.mo<-bj.plot[bj.plot$Locus %in% bj.mo$Locus,]
bj.mo.region<-data.frame(bj.mo$Chrom, as.numeric(bj.mo$Pos)-2500,
	as.numeric(bj.mo$Pos)+2500)
write.table(bj.mo.region,"../biallelic_outliers/rad_region/bjmo_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')
aj.bj<-aj.plot[aj.plot$Locus %in% aj.bj$Locus,]
aj.bj.region<-data.frame(aj.bj$Chrom, as.numeric(aj.bj$Pos)-2500,
	as.numeric(aj.bj$Pos)+2500)
write.table(aj.bj.region,"../biallelic_outliers/rad_region/ajbj_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')

sharedregion<-aj.plot[aj.plot$LocID %in% shared$LocID,]
shared.region<-data.frame(sharedregion$Chrom,as.numeric(sharedregion$Pos-2500),
	as.numeric(sharedregion$Pos)+2500)
write.table(shared.region,"../biallelic_outliers/rad_region/shared_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')

shared.chrom<-levels(as.factor(c(as.character(aj.out$Chrom),
	as.character(mo.out$Chrom),as.character(bj.out$Chrom))))
write.table(shared.chrom,"../biallelic_outliers/rad_region/top1_scaffolds.txt",
	quote=F,col.names=F,row.names=F,eol='\n')

###########LOOK INTO THE EXTREME OUTLIERS
#fm.extreme<-fm[fm$FEM.MAL >=0.11,]
mo.extreme<-mo[mo$FEM.MOM >= 0.08,c("Locus","Chrom","Pos","LocID","FEM.MOM")]
bj.extreme<-bj[bj$JUVIE.BREEDER >= 0.015,c("Locus","Chrom","Pos","LocID","JUVIE.BREEDER")]
aj.extreme<-aj[aj$ADULT.JUVIE >= 0.02,c("Locus","Chrom","Pos","LocID","ADULT.JUVIE")]

mo.ex.sum<-gw.sum[gw.sum$Locus %in% mo.extreme$Locus & 
	gw.sum$Pop %in% c("FEM","MOM"),]

bj.ex.sum<-gw.sum[gw.sum$Locus %in% bj.extreme$Locus & 
	gw.sum$Pop %in% c("JUVIE","BREEDER"),]

aj.ex.sum<-gw.sum[gw.sum$Locus %in% aj.extreme$Locus &
	gw.sum$Pop %in% c("ADULT","JUVIE"),]

par(mfrow=c(3,2))
hist(aj.ex.sum$Allele1Freq[aj.ex.sum$Pop=="ADULT"],xlab="",	ylab="",main="Adults")
hist(aj.ex.sum$Allele1Freq[aj.ex.sum$Pop=="JUVIE"],xlab="",	ylab="",main="Offspring")
hist(mo.extreme.sum$Allele1Freq[mo.extreme.sum$Pop=="MOM"],xlab="",	ylab="",main="Mothers")
hist(mo.extreme.sum$Allele1Freq[mo.extreme.sum$Pop=="FEM"],xlab="",	ylab="",main="Females")
hist(bj.ex.sum$Allele1Freq[bj.ex.sum$Pop=="JUVIE"],xlab="",	ylab="",main="Offspring")
hist(bj.ex.sum$Allele1Freq[bj.ex.sum$Pop=="BREEDER"],xlab="",	ylab="",main="Breeders")

snps<-read.table("../stacks/batch_1.catalog.snps.tsv",sep='\t',comment.char="#")
colnames(snps)<-c("SqlID","SampleID","LocusID","Column","Type","LR","Rank1",
	"Rank2","Rank3","Rank4")
tags<-read.table("../stacks/batch_1.catalog.tags.tsv",sep='\t',comment.char="#")
colnames(tags)<-c("SqlID","SampleID","LocusID","Chrom","BP","Strand",
	"SeqType","StackComponent","SeqID","Seq","Deleveraged","Blacklisted",
	"Lumberjackstack","LogLike")

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

col.order<-c("LocusID","Chrom","BP","NumSNPs","NumSNPsOut","Pop1N","Pop2N",
	"Pop1AF","Pop2AF","Pop1Ho","Pop2Ho","Pop1Hs","Pop2Hs","SNP1","SNP2",
	"SNP3","SNP4","Seq")

aj.ex.dat<-extract.info(aj.ex.sum,tags,snps)
aj.ex.dat<-aj.ex.dat[,col.order]
aj.ex.dat$Comparison<-"AO"
mo.ex.dat<-extract.info(mo.ex.sum,tags,snps)
mo.ex.dat<-mo.ex.dat[,col.order]
mo.ex.dat$Comparison<-"MF"
bj.ex.dat<-extract.info(bj.ex.sum,tags,snps)
bj.ex.dat<-bj.ex.dat[,col.order]
bj.ex.dat$Comparison<-"BO"

extreme.outliers<-rbind(aj.ex.dat,mo.ex.dat,bj.ex.dat)
write.csv(extreme.outliers,"ExtremeOutliers.csv",row.names=F)

aj.null<-aj.prune[!(aj.prune$Locus %in% aj.out$Locus),]
aj.null.sum<-gw.sum[gw.sum$Locus %in% aj.null$Locus,]
aj.null.sum<-aj.null.sum[aj.null.sum$Pop %in% c("ADULT","JUVIE"),]
aj.all<-extract.info(aj.null.sum,tags,snps)

bj.null<-bj.prune[!(bj.prune$Locus %in% bj.out$Locus),]
bj.null.sum<-gw.sum[gw.sum$Locus %in% bj.null$Locus,]
bj.null.sum<-bj.null.sum[bj.null.sum$Pop %in% c("BREEDER","JUVIE"),]
bj.all<-extract.info(bj.null.sum,tags,snps)

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
bj.out.dat<-gw.sum[gw.sum$Locus %in% bj.out$Locus,]
bj.out.dat<-bj.out.dat[bj.out.dat$Pop %in% c("BREEDER","JUVIE"),]
bj.out.dat<-extract.info(bj.out.dat,tags,snps)
nrad<-c(nrow(aj.null),nrow(aj.out.dat),nrow(aj.ex.dat),
	nrow(mo.all),nrow(mo.out.dat),nrow(mo.ex.dat),
	nrow(bj.all),nrow(bj.out.dat),nrow(bj.ex.dat))
library(scales)
png("SNPsPerRADLocus.png", height=7,width=7,units="in",res=300)
boxplot(aj.all$NumSNPs,aj.out.dat$NumSNPs,aj.ex.dat$NumSNPs,
	mo.all$NumSNPs,mo.out.dat$NumSNPs,mo.ex.dat$NumSNPs,
	bj.all$NumSNPs,bj.out.dat$NumSNPs,bj.ex.dat$NumSNPs,
	col=c("grey",alpha("green4",0.5),"green4",
		"grey",alpha("purple3",0.5),"purple3",
		"grey",alpha("dodgerblue",0.5),"dodgerblue"),
	names=rep("",9),ylim=c(0,70),axes=F)
axis(2,at=seq(0,70,10),las=1)
mtext("Number of SNPs Per RAD tag",2,line=2)
text(1:9,par("usr")[3] - 1, srt = 45, adj = 1,xpd = TRUE,
     labels = c("AO Null","AO Outliers","AO Extreme", "MF Null","MF Outliers",
		"MF Extreme","BO Null","BO Outliers","BO Extreme"))
text(1:9,y=66,labels=nrad)
dev.off()

snps.per.rad<-data.frame(rbind(cbind(aj.all$NumSNPs,"AO","Null"),
	cbind(aj.out.dat$NumSNPs,"AO","Outliers"),
	cbind(aj.ex.dat$NumSNPs,"AO","Extreme"),cbind(mo.all$NumSNPs,"FM","Null"),
	cbind(mo.out.dat$NumSNPs,"FM","Outliers"),
	cbind(mo.ex.dat$NumSNPs,"FM","Extreme"),cbind(bj.all$NumSNPs,"BO","Null"),
	cbind(bj.out.dat$NumSNPs,"BO","Outliers"),cbind(bj.ex.dat$NumSNPs,
	"BO","Extreme")))
colnames(snps.per.rad)<-c("NumSNPs","Comparison","SNPType")
spr.lm<-lm(as.numeric(NumSNPs)~Comparison+SNPType,data=snps.per.rad)
#> anova(spr.lm)
#Analysis of Variance Table
#
#Response: as.numeric(NumSNPs)
#              Df   Sum Sq Mean Sq F value    Pr(>F)    
#Comparison     2      877   438.5  2.2405    0.1064    
#SNPType        2    10365  5182.6 26.4778 3.206e-12 ***
#Residuals  59097 11567337   195.7     

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
	apply(bj.all[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(bj.out.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length)),
	do.call(cbind,lapply(
	apply(bj.ex.dat[,c("SNP1","SNP2","SNP3","SNP4")],2,grep,pattern="N"),
	length))))
Ns$Comparison<-c(rep("AO",3),rep("FM",3),rep("BO",3))
Ns$SNPType<-rep(c("Null","Outliers","Extreme"),3)
Ns$NumLoci<-c(nrow(aj.all),nrow(aj.out.dat),nrow(aj.ex.dat),
	nrow(mo.all),nrow(mo.out.dat),nrow(mo.ex.dat),
	nrow(bj.all),nrow(bj.out.dat),nrow(bj.ex.dat))
Ns$PropN<-rowSums(Ns[,1:4])/Ns$NumLoci

png("PropN.png",height=7,width=7,units="in",res=300)
bp<-barplot(Ns$PropN,
	col=c("grey",alpha("green4",0.5),"green4",
		"grey",alpha("purple3",0.5),"purple3",
		"grey",alpha("dodgerblue",0.5),"dodgerblue"),
	names="",ylim=c(0,0.5),las=1,ylab="Number of SNPs Per RAD tag")
text(bp,par("usr")[3] - 0.01, srt = 45, adj = 1,xpd = TRUE,
     labels = c("AO Null","AO Outliers","AO Extreme", "MF Null","MF Outliers",
		"MF Extreme","BO Null","BO Outliers","BO Extreme"))
text(bp,y=0.05,labels=Ns$NumLoci,srt=90)
dev.off()

#For outliers on loci with a bunch of SNPs, are most of the SNPs significant?
hist(aj.out.dat$NumSNPsOut/aj.out.dat$NumSNPs)
hist(bj.out.dat$NumSNPsOut/bj.out.dat$NumSNPs)
hist(mo.out.dat$NumSNPsOut/mo.out.dat$NumSNPs)


################Blast2Go Annotations
shared<-rbind(pj.fm,pj.mo,fm.mo)
shared<-shared[,c("Chrom","Pos","LocID","Locus")]
shared$comparison<-c(rep("PJ-FM",nrow(pj.fm)),rep("PJ-MO",nrow(pj.mo)),
	rep("FM-MO",nrow(fm.mo)))
shared2<-shared[shared$Chrom!="scaffold_985",]
shared2$start<-shared2$Pos-2500
shared2$end<-shared2$Pos+2500
shared2$start[shared2$start<0]<-0
shared2$SeqName<-paste(shared2$Chrom,"_",shared2$start,"-",shared2$end,sep="")
sharedb2g<-read.delim("rad_region/shared_region.blast2go.txt",
	header=T,sep='\t')
sh2<-merge(shared2,sharedb2g,by="SeqName",all=T)
write.table(sh2,"sharedintwo_blast2go.txt",row.names=F,col.names=T,quote=F,sep='\t')

#shared in three
sharedb2g[grep("scaffold_985",sharedb2g$SeqName),]#no hits

pj.unique<-pj.plot[pj.plot$Locus %in% pj.unique$Locus,
	c("Chrom","Pos","LocID","Locus")]
pj.unique$start<-pj.unique$Pos-2500
pj.unique$end<-pj.unique$Pos+2500
pj.unique$start[pj.unique$start < 0]<-"0"
pj.unique$SeqName<-paste(pj.unique$Chrom,"_",pj.unique$start,"-",
	pj.unique$end,sep="")
pjb2g<-read.delim("rad_region/pj_region.blast2go.txt",header=T,sep='\t')
pj.b2g<-merge(pj.unique,pjb2g,by="SeqName",all=T)
write.table(pj.b2g,"pj_region_withloc.blast2go.txt",row.names=F,col.names=T,sep='\t')

mo.unique<-mo.plot[mo.plot$Locus %in% mo.unique$Locus,
	c("Chrom","Pos","LocID","Locus")]
mo.unique$start<-mo.unique$Pos-2500
mo.unique$end<-mo.unique$Pos+2500
mo.unique$start[mo.unique$start < 0]<-"0"
mo.unique$SeqName<-paste(mo.unique$Chrom,"_",mo.unique$start,"-",
	mo.unique$end,sep="")
mob2g<-read.delim("rad_region/mo_region.blast2go.txt",header=T,sep='\t')
mo.b2g<-merge(mo.unique,mob2g,by="SeqName",all=T)
write.table(mo.b2g,"mo_region_withloc.blast2go.txt",row.names=F,col.names=T,sep='\t')



###########COMPARE TO PSTFST SIGNIFICANT LOCI
pstfst.bands<-scan(what="numeric",
	"E:/ubuntushare/popgen/sw_results/pstfst/sig_regions/Bands_extract.sh")
length(pstfst.scaff$V1[pstfst.scaff$V1 %in% fm.out$Chrom])
length(pstfst.scaff$V1[pstfst.scaff$V1 %in% mo.out$Chrom])
length(pstfst.scaff$V1[pstfst.scaff$V1 %in% pj.out$Chrom])
############################################################################
aj.ci<-c(mean(aj.prune$ADULT.JUVIE)+2.57583*sd(aj.prune$ADULT.JUVIE),
	mean(aj.prune$ADULT.JUVIE)-2.57583*sd(aj.prune$ADULT.JUVIE))
fm.ci<-c(mean(fm.prune$FEM.MAL)+2.57583*sd(fm.prune$FEM.MAL),
	mean(fm.prune$FEM.MAL)-(2.57583*sd(fm.prune$FEM.MAL)))
mo.ci<-c(mean(mo.prune$FEM.MOM)+(2.57583*sd(mo.prune$FEM.MOM)),
	mean(mo.prune$FEM.MOM)-(2.57583*sd(mo.prune$FEM.MOM)))

#top 5%
aj.plot<-aj.prune[order(aj.prune$ADULT.JUVIE),] #ascending
aj.top5<-c(aj.plot[round(nrow(aj.plot)*0.975),"ADULT.JUVIE"],
	aj.plot[round(nrow(aj.plot)*0.025),"ADULT.JUVIE"])
fm.plot<-fm.prune[order(fm.prune$FEM.MAL),]#ascending
fm.top5<-c(fm.plot[round(nrow(fm.plot)*0.975),"FEM.MAL"],
	fm.plot[round(nrow(fm.plot)*0.025),"FEM.MAL"])
mo.plot<-mo.prune[order(mo.prune$FEM.MOM),]#ascending
mo.top5<-c(mo.plot[round(nrow(mo.plot)*0.975),"FEM.MOM"],
	mo.plot[round(nrow(mo.plot)*0.025),"FEM.MOM"])
pj.plot<-pj.prune[order(pj.prune$JUVIE.PREGGER),]#ascending
pj.top5<-c(pj.plot[round(nrow(pj.plot)*0.975),"JUVIE.PREGGER"],
	pj.plot[round(nrow(pj.plot)*0.025),"JUVIE.PREGGER"])


#get model data
model<-read.delim("../sca_simulation_output/ddraddist.ss0.2alleles.error1.fst_out.txt")
model.pj<-model[model$POFst>0 & model$MaleAF < 0.95 & model$MaleAF > 0.05,]
pj.null<-c(mean(model.pj$POFst)+2.57583*sd(model.pj$POFst),
	mean(model.pj$POFst)-2.57583*sd(model.pj$POFst))

model.mo<-model[model$MFFst>0 & model$FemAF < 0.95 & model$FemAF > 0.05,]
mo.null<-c(mean(model.mo$MFFst)+2.57583*sd(model.mo$MFFst),
	mean(model.mo$MFFst)-(2.57583*sd(model.mo$MFFst)))

model.mf<-model[model$MFFst>0 & model$MaleAF < 0.95 & model$MaleAF > 0.05,]
mf.null<-c(mean(model.mf$MFFst)+(2.57583*sd(model.mf$MFFst)),
	mean(model.mf$MFFst)-(2.57583*sd(model.mf$MFFst)))

#plot
png("fst.biallelic.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(aj.prune, ci.dat=aj.ci,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.prune, ci.dat=fm.ci,fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.prune, ci.dat=mo.ci,fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

#plot with the model CIs
png("fst.biallelic.pruned.model.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(pj.prune, ci.dat=pj.null,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.prune, ci.dat=mf.null,fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.prune, ci.dat=mo.null,fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()



plot.fsts($Pos, aj$ADULT.JUVIE,pch=19,col="light grey")
points(aj.out$Pos,aj.out$ADULT.JUVIE,pch=19,col="dark grey")
points(aj.unique$Pos,aj.unique$ADULT.JUVIE,pch=19,col="dark green")

plot(fm$Pos, fm$FEM.MAL,pch=19,col="light grey")
points(fm.out$Pos,fm.out$FEM.MAL,pch=19,col="dark grey")
points(fm.unique$Pos,fm.unique$FEM.MAL,pch=19,col="dark green")

plot(mo$Pos, mo$FEM.MOM,pch=19,col="light grey")
points(mo.out$Pos,mo.out$FEM.MOM,pch=19,col="dark grey")
points(mo.unique$Pos,mo.unique$FEM.MOM,pch=19,col="dark green")

###top 1%
png("fst.biallelic.pruned.top1.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
aj<-plot.fsts(aj.prune, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
fm<-plot.fsts(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", 
	chrom.name="Chrom", axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
mo<-plot.fsts(mo.plot, ci.dat=c(mo.top1,0),fst.name="FEM.MOM", 
	chrom.name="Chrom", axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

aj.plot<-aj[order(aj$ADULT.JUVIE),] #ascending
aj.top1<-aj.plot[round(nrow(aj.plot)*0.99),"ADULT.JUVIE"]
aj.out1<-aj[aj$ADULT.JUVIE >= aj.top1,]
fm.plot<-fm[order(fm$FEM.MAL),]#ascending
fm.top1<-fm.plot[round(nrow(fm.plot)*0.99),"FEM.MAL"]
fm.out1<-fm[fm$FEM.MAL >= fm.top1,]
mo.plot<-mo[order(mo$FEM.MOM),]#ascending
mo.top1<-mo.plot[round(nrow(mo.plot)*0.99),"FEM.MOM"]
mo.out1<-mo[mo$FEM.MOM >= mo.top1,]

aj.un1<-aj.out1[!(aj.out1$Locus %in% fm.out1$Locus) &
	 !(aj.out1$Locus %in% mo.out1$Locus),]
fm.un1<-fm.out1[!(fm.out1$Locus %in% aj.out1$Locus) &
	!(fm.out1$Locus %in% mo.out1$Locus),]
mo.un1<-mo.out1[!(mo.out1$Locus %in% aj.out1$Locus) &
	!(mo.out1$Locus %in% fm.out1$Locus),]


png("biallelic.top1.png",width=7.5,height=10,units="in",res=300)
par(mfrow=c(3,1),lwd=1.3,las=1)
plot(aj$Pos, aj$ADULT.JUVIE,pch=19,col="light grey",axes=F,xlab="",
	ylab="")
points(aj.out1$Pos,aj.out1$ADULT.JUVIE,pch=19,col="dark blue")
points(aj.un1$Pos,aj.un1$ADULT.JUVIE,pch=19,col="dark green")
axis(2,pos=0)
mtext("Adults-Juveniles",3)
legend("topright",c("Shared","Unique"),col=c("dark blue","dark green"),
	pch=19,bty='n')
plot(fm$Pos, fm$FEM.MAL,pch=19,col="light grey",axes=F,xlab="",
	ylab="")
points(fm.out1$Pos,fm.out1$FEM.MAL,pch=19,col="dark blue")
points(fm.un1$Pos,fm.un1$FEM.MAL,pch=19,col="dark green")
axis(2,pos=0)
mtext("Females-Males",3)
mtext("Fst",2,las=0,outer=T)
plot(mo$Pos, mo$FEM.MOM,pch=19,col="light grey",axes=F,xlab="",
	ylab="")
points(mo.out1$Pos,mo.out1$FEM.MOM,pch=19,col="dark blue")
points(mo.un1$Pos,mo.un1$FEM.MOM,pch=19,col="dark green")
axis(2,pos=0)
mtext("Females-Mothers",3)
mtext("Location on genome",1,outer=T)
dev.off()

write.table(rownames(aj.un1),"unique.top1.aj.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(fm.un1),"unique.top1.fm.txt",quote=F,row.names=F,col.names=F)
write.table(rownames(mo.un1),"unique.top1.mo.txt",quote=F,row.names=F,col.names=F)

loci<-c(aj.out1$Locus,fm.out1$Locus,mo.out1$Locus)
locus.info<-strsplit(loci,split="\\.")
locus.info<-do.call("rbind",locus.info)
locus.info<-data.frame(cbind(loci,locus.info))
colnames(locus.info)<-c("Locus","Chrom","BP")
map<-read.table("../stacks/batch_1.plink.map")
map$Locus<-paste(map$V1,map$V4,sep=".")
stats<-read.table("../stacks/batch_1.sumstats.tsv")
stats$Locus<-paste(stats$V3,stats$V4,sep=".")

cat.loc<-stats[stats$Locus %in% locus.info$Locus,c("V2","Locus")]
cat.loc<-cat.loc[!duplicated(cat.loc$Locus),]
cat.loc<-merge(cat.loc,locus.info,by="Locus")
rad.loc<-cat.loc[!duplicated(cat.loc$V2),"V2"]
write.table(rad.loc,"top1.out.radloc.txt",quote=F,col.names=F,row.names=F)
write.table(levels(cat.loc$Chrom),"top1.scaffolds.txt",quote=F,col.names=F,row.names=F)
cat.loc$start<-as.numeric(as.character(cat.loc$BP))-2500
cat.loc$stop<-as.numeric(as.character(cat.loc$BP))+2500
cat.loc[cat.loc$start<0,"start"]<-0
write.table(cat.loc[,c("Chrom","start","stop")],"top1.2500bp.txt",
	quote=F,col.names=F,row.names=F)


