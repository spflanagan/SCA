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

##prune based on representation in the groups
sum.list<-split(gw.sum, gw.sum$Pop)
adt.n<-sum.list$ADULT[sum.list$ADULT$N > 440 & !is.na(sum.list$ADULT$Hs),]
juv.n<-sum.list$JUVIE[sum.list$JUVIE$N>200& !is.na(sum.list$JUVIE$Hs),]
fem.n<-sum.list$FEM[sum.list$FEM$N>130& !is.na(sum.list$FEM$Hs),]#130
mal.n<-sum.list$MAL[sum.list$MAL$N>200& !is.na(sum.list$MAL$Hs),]
mom.n<-sum.list$MOM[sum.list$MOM$N>100& !is.na(sum.list$MOM$Hs),]
pop.n<-sum.list$POP[sum.list$POP$N>87& !is.na(sum.list$POP$Hs),]
prg.n<-sum.list$PREGGER[sum.list$PREGGER$N>200& !is.na(sum.list$PREGGER$Hs),]
non.n<-sum.list$NONPREG[sum.list$NONPREG$N>14& !is.na(sum.list$NONPREG$Hs),]

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]
pop.n<-pop.n[pop.n$Allele1Freq > 0.05 & pop.n$Allele1Freq < 0.95,]
prg.n<-prg.n[prg.n$Allele1Freq > 0.05 & prg.n$Allele1Freq < 0.95,]
non.n<-non.n[non.n$Allele1Freq > 0.05 & non.n$Allele1Freq < 0.95,]


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

#get model data
model<-read.delim("../sca_simulation_output/knowndist.ss0.2alleles.fst_out.txt")
model<-read.delim("../sca_simulation_output/knowndist.ss0.2alleles.error01.fst_out.txt")
model.aj<-model[model$AOFst>0 & model$MaleAF < 0.95 & model$MaleAF > 0.05,]
aj.null<-c(mean(model.aj$AOFst)+2.57583*sd(model.aj$AOFst),
	mean(model.aj$AOFst)-2.57583*sd(model.aj$AOFst))

model.mo<-model[model$MDFst>0 & model$FemAF < 0.95 & model$FemAF > 0.05,]
mo.null<-c(mean(model.mo$MDFst)+2.57583*sd(model.mo$MDFst),
	mean(model.mo$MDFst)-(2.57583*sd(model.mo$MDFst)))

model.mf<-model[model$MFFst>0 & model$MaleAF < 0.95 & model$MaleAF > 0.05,]
mf.null<-c(mean(model.mf$MFFst)+(2.57583*sd(model.mf$MFFst)),
	mean(model.mf$MFFst)-(2.57583*sd(model.mf$MFFst)))

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
plot.fsts(fm.prune, ci.dat=mf.null,fst.name="MAL.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.prune, ci.dat=mo.null,fst.name="MOM.FEM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

##TROUBLESHOOTING MALE-FEMALE COMPARISON
 weirdos<-fm.prune[fm.prune$MAL.FEM > 0.03,]
weirdos<-weirdos[,c("Locus","Chrom","Pos","MAL.FEM")]
weirdsum<-gw.sum[gw.sum$Locus %in% weirdos$Locus,]
weirdsum<-weirdsum[weirdsum$Pop == "MAL" | weirdsum$Pop == "FEM",]

regsum<-gw.sum[gw.sum$Locus %in% fm.prune$Locus,] 
regsum<-regsum[regsum$Pop=="MAL" | regsum$Pop == "FEM",]
regsum<-regsum[!(regsum$Locus %in% weirdsum$Locus),]

png("HsInWeirdMalFemLoci.png",height=7,width=7,units="in",res=300)
par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
hist(weirdsum[weirdsum$Pop=="MAL","Hs"], breaks=50, main="",ylab="",xlab="")
legend("topright",bty="n","Weird Loci, Males")
hist(regsum[regsum$Pop=="MAL","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty="n","Regular Loci, Males")
hist(weirdsum[weirdsum$Pop=="FEM","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty='n',"Weird Loci, Females")
hist(regsum[regsum$Pop=="FEM","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty='n',"Regular Loci, Females")
mtext("Hs",1,outer=T)
mtext("Frequency",2,outer=T)
dev.off()


png("HvNInWeirdMalFemLoci.png",height=7,width=7,units="in",res=300)
par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(weirdsum[weirdsum$Pop=="MAL",c("N","Hs")],ylab="",xlab="",pch=15)
legend("bottomleft",bty="n","Hs Males",text.col="red")
plot(weirdsum[weirdsum$Pop=="MAL",c("N","Ho")],ylab="",xlab="",pch=15,col="blue")
legend("bottomleft",bty="n","Ho, Males",text.col="red")
plot(weirdsum[weirdsum$Pop=="FEM",c("N","Hs")],ylab="",xlab="",pch=19)
legend("bottomright",bty='n',"Hs, Females",text.col="red")
plot(weirdsum[weirdsum$Pop=="FEM",c("N","Ho")],ylab="",xlab="",pch=19,col="blue")
legend("bottomleft",bty='n',"Ho, Females",text.col="red")
mtext("Heterozygosity",2,outer=T)
mtext("N",1,outer=T)
mtext("Weird Loci",3,outer=T)
dev.off()

 


fem.n<-sum.list$FEM[!is.na(sum.list$FEM$Hs),]
 plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#OK
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#OK
fem.n<-fem.n[fem.n$N>100,]
plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#goodbye

removedfem<-fem.n[fem.n$N < 100,]
plot(gw.fst[gw.fst$Locus %in% removedfem$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])
rem.fem<-gw.fst[gw.fst$Locus %in% removedfem$Locus & gw.fst$MAL.FEM > 0.025 
	& gw.fst$MAL.FEM < 0.05,c("Locus","Chrom","Pos","MAL.FEM")]
rem.fem<-merge(rem.fem, removedfem,by="Locus")

