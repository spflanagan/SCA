#Author: Sarah P. Flanagan
#Last Updated: 27 April 2016
#Started Date: 11 February 2016
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

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]

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
mo.prune<-mo.prune[mo.prune$FEM.MOM>0,]


write.table(pj.prune, "pj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(fm.prune, "fm.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(mo.prune, "mo.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(aj.prune,"aj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')
write.table(bj.prune,"bj.plot.txt",row.names=F,col.names=T,quote=F,sep='\t')

aj.plot<-aj.prune[order(aj.prune$ADULT.JUVIE),] #ascending
aj.top1<-aj.plot[round(nrow(aj.plot)*0.99),"ADULT.JUVIE"]
aj.out1<-aj.plot[aj.plot$ADULT.JUVIE >= aj.top1,]
fm.plot<-fm.prune[order(fm.prune$FEM.MAL),]#ascending
fm.top1<-fm.plot[round(nrow(fm.plot)*0.99),"FEM.MAL"]
fm.out1<-fm.prune[fm.prune$FEM.MAL >= fm.top1,]
mo.plot<-mo.prune[order(mo.prune$FEM.MOM),]#ascending
mo.top1<-mo.plot[round(nrow(mo.plot)*0.99),"FEM.MOM"]
mo.out1<-mo.prune[mo.prune$FEM.MOM >= mo.top1,]
#pj.plot<-pj.prune[order(pj.prune$JUVIE.PREGGER),] #ascending
#pj.top1<-pj.plot[round(nrow(pj.plot)*0.99),"JUVIE.PREGGER"]
#pj.out1<-pj.prune[pj.prune$JUVIE.PREGGER >= pj.top1,]
#bj.plot<-bj.prune[order(bj.prune$JUVIE.BREEDER),]#ascending
#bj.top1<-bj.plot[round(nrow(bj.plot)*0.99),"JUVIE.BREEDER"]
#bj.out1<-bj.prune[bj.prune$JUVIE.BREEDER >= bj.top1,]

#plot with the top1%
#png("fst.top1.comp3_redo.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)

mo<-plot.fsts(mo.plot, ci.dat=c(mo.top1,0),fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)

fm<-plot.fsts(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
legend("top","Male-Female", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)

aj<-plot.fsts(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", 
	chrom.name="Chrom", axis.size=0.75,bp.name="Pos",sig.col=c("red","black"))
legend("top","Adult-Offspring", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)


#dev.off()



###COMPARISONS
aj.out<-aj[aj$ADULT.JUVIE >= aj.top1[1],]
fm.out<-fm[fm$FEM.MAL >= fm.top1[1],]
mo.out<-mo[mo$FEM.MOM >= mo.top1[1],]
#pj.out<-pj[pj$JUVIE.PREGGER >= pj.top1[1],]
#bj.out<-bj[bj$JUVIE.BREEDER >= bj.top1[1],]

aj.unique<-aj.out[!(aj.out$Locus %in% fm.out$Locus) & 
	!(aj.out$Locus %in% mo.out$Locus),]
fm.unique<-fm.out[!(fm.out$Locus %in% aj.out$Locus) & 
	!(fm.out$Locus %in% mo.out$Locus),]
mo.unique<-mo.out[!(mo.out$Locus %in% aj.out$Locus) &
	!(mo.out$Locus %in% fm.out$Locus),]
#pj.unique<-pj.out[#!(pj.out$Locus %in% aj.out$Locus) &
#	!(pj.out$Locus %in% fm.out$Locus) & !(pj.out$Locus %in% mo.out$Locus),]
bj.unique<-bj.out[!(bj.out$Locus %in% aj.out$Locus) &
	!(bj.out$Locus %in% mo.out$Locus),]
shared<-aj.out[(aj.out$LocID %in% mo.out$LocID) & 
	(aj.out$LocID %in% fm.out$LocID),]

aj.fm<-aj.out[(aj.out$LocID %in% fm.out$LocID),]
aj.mo<-aj.out[(aj.out$LocID %in% mo.out$LocID),]
fm.mo<-fm.out[fm.out$LocID %in% mo.out$LocID,]
length(levels(factor(c(as.character(aj.fm$Locus),as.character(aj.mo$Locus),
	as.character(fm.mo$Locus)))))

png("fst.selection.episodes_redo.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
mo<-plot.fsts(mo.plot, ci.dat=c(mo.top1,0),fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0,bp.name="Pos",sig.col=c("purple3","black"))#ignore error
points(mo$Pos[mo$LocID %in% shared$LocID & mo$FEM.MOM >= mo.top1],
	mo$FEM.MOM[mo$LocID %in% shared$LocID& mo$FEM.MOM >= mo.top1],
	col="red",pch=8)
axis(2,at=seq(0,0.15,0.05),pos=0,las=1,cex.axis=0.75)
legend("top","Females-Inferred Mothers",bty='n',cex=0.75,text.font=2)

fm<-plot.fsts(fm.plot, ci.dat=c(fm.top1,0),fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("green4","black"))
points(fm$Pos[fm$LocID %in% shared$LocID & fm$FEM.MAL >= fm.top1],
	fm$FEM.MAL[fm$LocID %in% shared$LocID& fm$FEM.MAL >= fm.top1],
	col="red",pch=8)
legend("top","Male-Female",bty='n',cex=0.75,text.font=2)
#mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1.5)

aj<-plot.fsts(aj.plot, ci.dat=c(aj.top1,0),fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0, bp.name="Pos",sig.col=c("dodgerblue","black"))#ignore error
axis(2,at=c(0,0.025,0.05),pos=0,las=1,cex.axis=0.75)
points(aj$Pos[aj$LocID %in% shared$LocID& aj$ADULT.JUVIE >= aj.top1],
	aj$ADULT.JUVIE[aj$LocID %in% shared$LocID& aj$ADULT.JUVIE >= aj.top1],
	col="red",pch=8)
legend("top","Adult-Offspring",bty='n',cex=0.75,text.font=2)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1,las=0)

par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","purple3","dodgerblue","red"),pch=c(19,19,19,8),
	c("Viability Selection","Sexual Selection","Overall Selection","Shared in all"),
	bg="white",ncol=4,box.lty=0)
dev.off()

##WRITE TO FILE
write.table(levels(as.factor(aj.unique$LocID)),
	"../biallelic_outliers/AJ_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(fm.unique$LocID)),
	"../biallelic_outliers/FM_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(mo.unique$LocID)),
	"../biallelic_outliers/MO_1outliers.txt",quote=F,
	col.names=F,row.names=F,sep='\t')
#write.table(levels(as.factor(pj.unique$LocID)),
#	"../biallelic_outliers/PJ_1outliers.txt",quote=F,
#	col.names=F,row.names=F,sep='\t')
#write.table(levels(as.factor(bj.out$LocID)),
#	"../biallelic_outliers/BJ_outliers.txt",quote=F,
#	col.names=F,row.names=F,sep='\t')
write.table(levels(as.factor(shared$LocID)),
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

aj.unique<-aj.plot[aj.plot$Locus %in% aj.unique$Locus,]
aj.rad.region<-data.frame(aj.unique$Chrom,as.numeric(aj.unique$Pos)-2500,
	as.numeric(aj.unique$Pos)+2500)
write.table(aj.rad.region,"../biallelic_outliers/rad_region/aj_extract.sh",
	quote=F,col.names=F,row.names=F,sep="\t")
fm.unique<-fm.plot[fm.plot$Locus %in% fm.unique$Locus,]
fm.rad.region<-data.frame(fm.unique$Chrom,as.numeric(fm.unique$Pos)-2500,
	as.numeric(fm.unique$Pos)+2500)
write.table(fm.rad.region,"../biallelic_outliers/rad_region/fm_1extract_redo.sh",
	quote=F,col.names=F,row.names=F,sep="\t")
mo.unique<-mo.plot[mo.plot$Locus %in% mo.unique$Locus,]
mo.rad.region<-data.frame(mo.unique$Chrom,as.numeric(mo.unique$Pos)-2500,
	as.numeric(mo.unique$Pos)+2500)
write.table(mo.rad.region,"../biallelic_outliers/rad_region/mo_1extract.sh",
	quote=F,col.names=F,row.names=F,sep="\t")
#bj.unique<-bj.plot[bj.plot$Locus %in% bj.unique$Locus,]
#bj.rad.region<-data.frame(bj.unique$Chrom,as.numeric(bj.unique$Pos)-2500,
#	as.numeric(bj.unique$Pos)+2500)
#write.table(bj.rad.region,"../biallelic_outliers/rad_region/bj_1extract.sh",
#	quote=F,col.names=F,row.names=F,sep="\t")

aj.mo<-aj.plot[aj.plot$Locus %in% aj.mo$Locus,]
aj.mo.region<-data.frame(aj.mo$Chrom, as.numeric(aj.mo$Pos)-2500,
	as.numeric(aj.mo$Pos)+2500)
write.table(aj.mo.region,"../biallelic_outliers/rad_region/ajmo_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')
fm.mo<-fm.plot[fm.plot$Locus %in% fm.mo$Locus,]
fm.mo.region<-data.frame(fm.mo$Chrom, as.numeric(fm.mo$Pos)-2500,
	as.numeric(fm.mo$Pos)+2500)
write.table(fm.mo.region,"../biallelic_outliers/rad_region/fmmo_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')
aj.fm<-aj.plot[aj.plot$Locus %in% aj.fm$Locus,]
aj.fm.region<-data.frame(aj.fm$Chrom, as.numeric(aj.fm$Pos)-2500,
	as.numeric(aj.fm$Pos)+2500)
write.table(aj.fm.region,"../biallelic_outliers/rad_region/ajfm_extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')

sharedregion<-aj.plot[aj.plot$LocID %in% shared$LocID,]
shared.region<-data.frame(sharedregion$Chrom,as.numeric(sharedregion$Pos-2500),
	as.numeric(sharedregion$Pos)+2500)
write.table(shared.region,"../biallelic_outliers/rad_region/shared_2extract.sh",
	quote=F,col.names=F,row.names=F,sep='\t')

shared.chrom<-levels(as.factor(c(as.character(aj.out$Chrom),
	as.character(mo.out$Chrom),as.character(bj.out$Chrom))))
write.table(shared.chrom,"../biallelic_outliers/rad_region/top1_scaffolds.txt",
	quote=F,col.names=F,row.names=F,eol='\n')

###########LOOK INTO THE EXTREME OUTLIERS
fm.extreme<-fm[fm$FEM.MAL >=0.11,c("Locus","Chrom","Pos","LocID","FEM.MAL")]
mo.extreme<-mo[mo$FEM.MOM >= 0.08,c("Locus","Chrom","Pos","LocID","FEM.MOM")]
#bj.extreme<-bj[bj$JUVIE.BREEDER >= 0.015,c("Locus","Chrom","Pos","LocID","JUVIE.BREEDER")]
aj.extreme<-aj[aj$ADULT.JUVIE >= 0.02,c("Locus","Chrom","Pos","LocID","ADULT.JUVIE")]

mo.ex.sum<-gw.sum[gw.sum$Locus %in% mo.extreme$Locus & 
	gw.sum$Pop %in% c("FEM","MOM"),]

fm.ex.sum<-gw.sum[gw.sum$Locus %in% fm.extreme$Locus & 
	gw.sum$Pop %in% c("FEM","MAL"),]

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
aj.ex.dat$Comparison<-"Adult-Offspring"
mo.ex.dat<-extract.info(mo.ex.sum,tags,snps)
mo.ex.dat<-mo.ex.dat[,col.order]
mo.ex.dat$Comparison<-"Mothers-Females"
fm.ex.dat<-extract.info(fm.ex.sum,tags,snps)
fm.ex.dat<-fm.ex.dat[,col.order]
fm.ex.dat$Comparison<-"Female-Male"

extreme.outliers<-rbind(aj.ex.dat,mo.ex.dat,fm.ex.dat)
write.csv(extreme.outliers,"ExtremeOutliers.csv",row.names=F)

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


################Blast2Go Annotations
setwd("../biallelic_outliers/blastresults")
out.cnames<-c("LocID","NumSNPs","Chrom.x","Pos","Description","Length",
	"X.Hits","e.Value","sim.mean","X.GO","GO.Names.list","Enzyme.Codes.list","Seq")
out.names<-c("LocID","NumSNPs","Chrom","Pos","Description","Length",
	"NumHits","e.Value","sim.mean","NumGO","GO.Names.list","Enzyme.Codes.list","Seq")
blast2go.files<-list.files(pattern="blast2go")
ajfm.b2g<-read.delim("ajbj_blast2go.txt",sep='\t',header=T)
ajfm.b2g$Comparison<-"AO-BO"
ajmo.b2g<-read.delim("ajmo_blast2go.txt",sep='\t',header=T)
ajmo.b2g$Comparison<-"AO-FM"
fmmo.b2g<-read.delim("bjmo_blast2go.txt",sep='\t',header=T)
fmmo.b2g$Comparison<-"BO-FM"
shared.b2g<-read.delim("shared_blast2go.txt",sep='\t',header=T)

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
#fmmo<-merge(fm.mo,bjmo.b2g,by="SeqName")
#fmmo<-merge(fmmo,mo.out.dat,by.x="LocID",by.y="LocusID")
#fmmo<-fmmo[,out.cnames]
#colnames(fmmo)<-out.names
#shared<-rbind(ajbj,ajmo,bjmo)

shared<-rbind(aj.fm,aj.mo,fm.mo)
shared.blast<-merge(shared,shared.b2g,by="SeqName")
all.out.dat<-rbind(mo.out.dat,aj.out.dat,fm.out.dat)
all.out.dat<-all.out.dat[!duplicated(all.out.dat$LocusID),]
shared.blast<-merge(shared.blast,all.out.dat,by.x="LocID",by.y="LocusID")
shared.blast<-cbind(shared.blast[,"Comparison"],shared.blast[,out.cnames])
colnames(shared.blast)<-c("Comparison",out.names)

write.table(shared.blast,"S1_SharedOutliersBlast.txt",sep='\t',col.names=T,
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
aj.blast<-merge(aj.unique,aj.b2g,by="SeqName",all=T)
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
fm.blast<-merge(fm.unique,fm.b2g,by="SeqName",all=T)
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
mo.b2g<-merge(mo.unique,mob2g,by="SeqName",all=T)
mo.blast<-merge(mo.b2g,mo.out.dat,by.x="LocID",by.y="LocusID")
mo.blast<-mo.blast[,out.cnames]
colnames(mo.blast)<-out.names
mo.blast$Comparison<-"Mothers-Females"

unique<-rbind(aj.blast,fm.blast,mo.blast)
write.table(unique,"S2_UniqueOutliersBlast.txt",row.names=F,col.names=T,
	sep='\t',quote=F)
###########COMPARE TO PSTFST SIGNIFICANT LOCI
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


