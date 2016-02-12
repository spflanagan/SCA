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
adt.grp.trim<-rownames(gw.loc.count[gw.loc.count$ADULTCount>440,])
juv.grp.trim<-rownames(gw.loc.count[gw.loc.count$JUVIECount>160,])
fem.grp.trim<-rownames(gw.loc.count[gw.loc.count$FEMCount>130,])
mal.grp.trim<-rownames(gw.loc.count[gw.loc.count$MALCount>196,])
mom.grp.trim<-rownames(gw.loc.count[gw.loc.count$MOMCount>10,])
pop.grp.trim<-rownames(gw.loc.count[gw.loc.count$POPCount>87,])
prg.grp.trim<-rownames(gw.loc.count[gw.loc.count$PREGGERCount>188,])
non.grp.trim<-rownames(gw.loc.count[gw.loc.count$NONPREGCount>8,])

#the comparisons
aj.loci<-rownames(gw.loc.count[gw.loc.count$ADULTCount>440 & 
	gw.loc.count$JUVIECount>160,])
fm.loci<-rownames(gw.loc.count[gw.loc.count$FEMCount>130 & 
	gw.loc.count$MALCount>196,])
mo.loci<-rownames(gw.loc.count[gw.loc.count$POPCount>87 & 
	gw.loc.count$MOMCount>10,])
np.loci<-rownames(gw.loc.count[gw.loc.count$PREGGERCount>188 & 
	gw.loc.count$NONPREGCount>8,])

#apply both filters
adt.prune<-adt.af.trim[names(adt.af.trim) %in% adt.grp.trim]
juv.prune<-juv.af.trim[names(juv.af.trim) %in% juv.grp.trim]
fem.prune<-fem.af.trim[names(fem.af.trim) %in% fem.grp.trim]
mal.prune<-mal.af.trim[names(mal.af.trim) %in% mal.grp.trim]
mom.prune<-mom.af.trim[names(mom.af.trim) %in% mom.grp.trim]
pop.prune<-pop.af.trim[names(pop.af.trim) %in% pop.grp.trim]
prg.prune<-prg.af.trim[names(prg.af.trim) %in% prg.grp.trim]
non.prune<-non.af.trim[names(non.af.trim) %in% non.grp.trim]
#comparisons
aj.prune<-aj.loci[aj.loci %in% aj.maf]
fm.prune<-fm.loci[fm.loci %in% fm.maf]
mo.prune<-mo.loci[mo.loci %in% mo.maf]
np.prune<-np.loci[np.loci %in% np.maf]

all.pruned<-gw.alleles[gw.alleles$Locus %in% aj.prune &
	gw.alleles$Locus %in% fm.prune &
	gw.alleles$Locus %in% mo.prune,]

gw.plot<-data.frame(Locus=gw.fst$Locus,Adult.Juvie=gw.fst$ADULT.JUVIE, 
	Fem.Mal=gw.fst$FEM.MAL, Fem.Mom=gw.fst$POP.MOM, 
	Nonpreg.Pregger=gw.fst$NONPREG.PREGGER)
gw.loc.info<-data.frame(Locus=gw.cat[,3],Chrom=gw.cat[,4],BP=gw.cat[,5])
gw.plot<-merge(gw.plot, gw.loc.info, by.x="Locus",by.y="Locus")

aj.plot<-gw.plot[gw.plot$Locus %in% aj.prune,]
write.table(aj.plot, "aj.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
fm.plot<-gw.plot[gw.plot$Locus %in% fm.prune,]
write.table(fm.plot, "fm.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
mo.plot<-gw.plot[gw.plot$Locus %in% mo.prune,]
write.table(gw.plot, "mo.plot.txt",row.names=F,col.names=F,quote=F,sep='\t')
np.plot<-gw.plot[gw.plot$Locus %in% np.prune,]



aj.ci<-c(mean(aj.plot$Adult.Juvie)+2.57583*sd(aj.plot$Adult.Juvie),
	mean(aj.plot$Adult.Juvie)-2.57583*sd(aj.plot$Adult.Juvie))
fm.ci<-c(mean(fm.plot$Fem.Mal)+2.57583*sd(fm.plot$Fem.Mal),
	mean(fm.plot$Fem.Mal)-(2.57583*sd(fm.plot$Fem.Mal)))
mo.ci<-c(mean(mo.plot$Fem.Mom)+(2.57583*sd(mo.plot$Fem.Mom)),
	mean(mo.plot$Fem.Mom)-(2.57583*sd(mo.plot$Fem.Mom)))
np.ci<-c(mean(np.plot$Nonpreg.Pregger)+(2.57583*sd(np.plot$Nonpreg.Pregger)),
	mean(np.plot$Nonpreg.Pregger)-(2.57583*sd(np.plot$Nonpreg.Pregger)))


