#Author: Sarah P. Flanagan
#Date: 11 February 2016
#Purpose: Analyze the output from gwsca_biallelic.

setwd("E:/ubuntushare/SCA/results/biallelic")
gw.fst<-read.delim("gwsca_fsts.txt")
gw.sum<-read.delim("gwsca_summary.txt")

gw.fst<-gw.fst[,c("Chrom","Pos","ADULT.JUVIE","MAL.FEM","PREGGER.OFF",
	"MOM.FEM")]
gw.fst$Locus<-paste(gw.fst$Chrom,gw.fst$Pos,sep=".")
gw.sum$Locus<-paste(gw.sum$Chrom,gw.sum$Pos,sep=".")


##NOW prune based on representation in the groups
sum.list<-split(gw.sum, gw.sum$Pop)
adt.n<-sum.list$ADULT[sum.list$ADULT$N > 440 & !is.na(sum.list$ADULT$Hs),]
juv.n<-sum.list$JUVIE[sum.list$JUVIE$N>160& !is.na(sum.list$JUVIE$Hs),]
fem.n<-sum.list$FEM[sum.list$FEM$N>130& !is.na(sum.list$FEM$Hs),]
mal.n<-sum.list$MAL[sum.list$MAL$N>196& !is.na(sum.list$MAL$Hs),]
mom.n<-sum.list$MOM[sum.list$MOM$N>100& !is.na(sum.list$MOM$Hs),]
pop.n<-sum.list$POP[sum.list$POP$N>87& !is.na(sum.list$POP$Hs),]
prg.n<-sum.list$PREGGER[sum.list$PREGGER$N>188& !is.na(sum.list$PREGGER$Hs),]
non.n<-sum.list$NONPREG[sum.list$NONPREG$N>8& !is.na(sum.list$NONPREG$Hs),]


#comparisons
aj.prune<-gw.fst[gw.fst$Locus %in% adt.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
aj.prune<-aj.prune[aj.prune$ADULT.JUVIE>0,]
fm.prune<-gw.fst[gw.fst$Locus %in% mal.n$Locus&gw.fst$Locus %in% fem.n$Locus, ]
fm.prune<-fm.prune[fm.prune$MAL.FEM>0,]
mo.prune<-gw.fst[gw.fst$Locus %in% fem.n$Locus&gw.fst$Locus %in% mom.n$Locus, ]
mo.prune<-mo.prune[mo.prune$MOM.FEM>0,]

write.table(aj.plot, "aj.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
write.table(fm.plot, "fm.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
write.table(gw.plot, "mo.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')


gw.plot<-gw.fst[gw.fst$Locus %in% aj.prune$Locus |
	gw.fst$Locus %in% fm.prune$Locus |
	gw.fst$Locus %in% mo.prune$Locus,]

aj.ci<-c(mean(aj.prune$ADULT.JUVIE)+2.57583*sd(aj.prune$ADULT.JUVIE),
	mean(aj.prune$ADULT.JUVIE)-2.57583*sd(aj.prune$ADULT.JUVIE))
fm.ci<-c(mean(fm.prune$MAL.FEM)+2.57583*sd(fm.prune$MAL.FEM),
	mean(fm.prune$MAL.FEM)-(2.57583*sd(fm.prune$MAL.FEM)))
mo.ci<-c(mean(mo.prune$MOM.FEM)+(2.57583*sd(mo.prune$MOM.FEM)),
	mean(mo.prune$MOM.FEM)-(2.57583*sd(mo.prune$MOM.FEM)))

model.aj<-read.delim("../sca_simulation_output/knowndist.ss0.2allelesao.txt")
model.aj<-model.aj[model.aj$Fst>0,]
aj.null<-c(mean(model.aj$Fst)+2.57583*sd(model.aj$Fst),
	mean(model.aj$Fst)-2.57583*sd(model.aj$Fst))

model.mo<-read.delim("../sca_simulation_output/knowndist.ss0.2allelesgp.txt")
model.mo<-model.mo[model.mo$Fst>0,]
fm.null<-c(mean(model.mo$Fst)+2.57583*sd(model.mo$Fst),
	mean(model.mo$Fst)-(2.57583*sd(model.mo$Fst)))

model.mf<-read.delim("../sca_simulation_output/knowndist.ss0.2allelesmf.txt")
model.mf<-model.mf[model.mf$Fst>0,]
mo.null<-c(mean(model.mf$Fst)+(2.57583*sd(model.mf$Fst)),
	mean(model.mf$Fst)-(2.57583*sd(model.mf$Fst)))

#plot
png("fst.biallelic.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(aj.prune, ci.dat=aj.ci,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.prune, ci.dat=fm.ci,fst.name="MAL.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.prune, ci.dat=mo.ci,fst.name="MOM.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

#plot with the model CIs
png("fst.biallelic.model.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(aj.prune, ci.dat=aj.null,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.prune, ci.dat=fm.null,fst.name="MAL.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.prune, ci.dat=mo.null,fst.name="MOM.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()
