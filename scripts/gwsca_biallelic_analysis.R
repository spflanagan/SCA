#Author: Sarah P. Flanagan
#Date: 11 February 2016
#Purpose: Analyze the output from gwsca_biallelic.

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")
gw.fst<-read.delim("gwsca_fsts.txt")
gw.sum<-read.delim("gwsca_summary.txt")
#gw.fst<-gw.fst[,c("Chrom","Pos","ADULT.JUVIE","MAL.FEM","PREGGER.OFF",
#	"MOM.FEM")]
gw.fst$Locus<-paste(gw.fst$Chrom,gw.fst$Pos,sep=".")
gw.sum$Locus<-paste(gw.sum$Chrom,gw.sum$Pos,sep=".")

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
#pop.n<-sum.list$POP[sum.list$POP$N>57& !is.na(sum.list$POP$Hs),]
prg.n<-sum.list$PREGGER[sum.list$PREGGER$N>159& !is.na(sum.list$PREGGER$Hs),]
non.n<-sum.list$NONPREG[sum.list$NONPREG$N>14& !is.na(sum.list$NONPREG$Hs),]

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]
#pop.n<-pop.n[pop.n$Allele1Freq > 0.05 & pop.n$Allele1Freq < 0.95,]
prg.n<-prg.n[prg.n$Allele1Freq > 0.05 & prg.n$Allele1Freq < 0.95,]
non.n<-non.n[non.n$Allele1Freq > 0.05 & non.n$Allele1Freq < 0.95,]


#comparisons
#viability
aj.prune<-gw.fst[gw.fst$Locus %in% adt.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
aj.prune<-aj.prune[aj.prune$ADULT.JUVIE>0,]
fm.prune<-gw.fst[gw.fst$Locus %in% mal.n$Locus&gw.fst$Locus %in% fem.n$Locus, ]
fm.prune<-fm.prune[fm.prune$FEM.MAL>0,]
#sexual
mo.prune<-gw.fst[gw.fst$Locus %in% fem.n$Locus&gw.fst$Locus %in% mom.n$Locus, ]
mo.prune<-mo.prune[mo.prune$FEM.MOM>0,]
#gametic
pj.prune<-gw.fst[gw.fst$Locus %in% prg.n$Locus&gw.fst$Locus %in% juv.n$Locus, ]
pj.prune<-pj.prune[pj.prune$JUVIE.PREGGER>0,]


write.table(aj.prune, "aj.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
write.table(fm.prune, "fm.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
write.table(mo.prune, "mo.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')


#gw.plot<-gw.fst[gw.fst$Locus %in% aj.prune$Locus &
#	gw.fst$Locus %in% fm.prune$Locus &
#	gw.fst$Locus %in% mo.prune$Locus,]

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
model<-read.delim("../sca_simulation_output/ddraddist.ss0.2alleles.fst_out.txt")
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

#plot with the top5%
png("fst.viability.top5.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
aj<-plot.fsts(aj.plot, ci.dat=aj.top5,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
fm<-plot.fsts(fm.plot, ci.dat=fm.top5,fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Male-Female", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

png("fst.sexual.top5.png",height=150,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
mo<-plot.fsts(mo.plot, ci.dat=mo.top5,fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

png("fst.gametic.top5.png",height=150,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
pj<-plot.fsts(pj.plot, ci.dat=pj.top5,fst.name="JUVIE.PREGGER", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()



###COMPARISONS
aj.out<-aj[aj$ADULT.JUVIE >= aj.top5[1] |	aj$ADULT.JUVIE <= aj.top5[2],]
fm.out<-fm[fm$FEM.MAL >= fm.top5[1] | fm$FEM.MAL <= fm.top5[2],]
mo.out<-mo[mo$FEM.MOM >= mo.top5[1] | mo$FEM.MOM <= mo.top5[2],]
pj.out<-pj[pj$JUVIE.PREGGER >= pj.top5[1] | pj$JUVIE.PREGGER <= pj.top5[2],]
vi.out<-aj.out[aj.out$Locus %in% fm.out$Locus,]

aj.unique<-aj.out[!(aj.out$Locus %in% fm.out$Locus) & 
	!(aj.out$Locus %in% mo.out$Locus) & !(aj.out$Locus %in% pj.out$Locus),]
fm.unique<-fm.out[!(fm.out$Locus %in% aj.out$Locus) & 
	!(fm.out$Locus %in% mo.out$Locus) & !(fm.out$Locus %in% pj.out$Locus),]
mo.unique<-mo.out[!(mo.out$Locus %in% aj.out$Locus) &
	!(mo.out$Locus %in% fm.out$Locus) & !(mo.out$Locus %in% pj.out$Locus),]
pj.unique<-pj.out[!(pj.out$Locus %in% aj.out$Locus) &
	!(pj.out$Locus %in% fm.out$Locus) & !(pj.out$Locus %in% mo.out$Locus),]

shared<-aj.out[(aj.out$Locus %in% fm.out$Locus) & 
	(aj.out$Locus %in% mo.out$Locus) & (aj.out$Locus %in% pj.out$Locus),]


png("fst.viability.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(2,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(aj.plot, ci.dat=aj.top5,fst.name="ADULT.JUVIE", chrom.name="Chrom"
	, axis.size=0.75, bp.name="Pos",sig.col=c("green4","green4"))
points(vi.out$Pos,vi.out$ADULT.JUVIE,col="dodgerblue",pch=16)
points(shared$Pos,shared$ADULT.JUVIE,col="red",pch=8)
legend("top","Adult-Juvenile",bty='n',cex=0.75)
plot.fsts(fm.plot, ci.dat=fm.top5,fst.name="FEM.MAL", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("green4","green4"))
points(vi.out$Pos,vi.out$FEM.MAL,col="dodgerblue",pch=16)
points(shared$Pos,shared$FEM.MAL,col="red",pch=8)
legend("top","Male-Female",bty='n',cex=0.75)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("green4","dodgerblue","red"),pch=c(19,16,8),
	c("5% Outliers","Viability Selection","Shared in all"),
	bg="white",ncol=3,box.lty=0)
dev.off()

png("fst.sexual.png",height=150,width=300,units="mm",res=300)
par(oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
mo<-plot.fsts(mo.plot, ci.dat=mo.top5,fst.name="FEM.MOM", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("purple3","purple3"))
points(shared$Pos,shared$FEM.MOM,col="red",pch=8)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("purple3","red"),pch=c(19,8),
	c("5% Outliers","Shared in all"),bg="white",ncol=2,box.lty=0)
dev.off()

png("fst.gametic.png",height=150,width=300,units="mm",res=300)
par(oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
pj<-plot.fsts(pj.plot, ci.dat=pj.top5,fst.name="JUVIE.PREGGER", chrom.name="Chrom"
	, axis.size=0.75,bp.name="Pos",sig.col=c("goldenrod","goldenrod"))
points(shared$Pos,shared$JUVIE.PREGGER,col="red",pch=8)
mtext("Genomic Location", 1, outer=T, cex=1)
mtext(expression(italic(F)[italic(ST)]), 2, outer=T, cex=1)
par(fig = c(0, 1, 0, 1), oma=c(2,1,0,1), mar = c(0, 0, 0, 0), new = TRUE,
	cex=1)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",col=c("goldenrod","red"),pch=c(19,8),
	c("5% Outliers","Shared in all"),bg="white",ncol=2,box.lty=0)
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


