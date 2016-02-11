#Author: Sarah P. Flanagan
#Date: 11 February 2016
#Purpose: Analyze the output from gwsca_biallelic.

setwd("E:/ubuntushare/SCA/results/biallelic")
gw.fst<-read.delim("gwsca_fsts.txt")
gw.sum<-read.delim("gwsca_summary.txt")


gw.all.loc<-split(gw.alleles, gw.alleles$Locus)
gw.all.sum<-lapply(gw.all.loc, function(x) colSums(x[,3:11]))
gw.loc.count<-data.frame(matrix(unlist(gw.all.sum), nrow=length(gw.all.sum), byrow=T),stringsAsFactors=FALSE)
colnames(gw.loc.count)<-names(gw.all.sum[[1]])
rownames(gw.loc.count)<-names(gw.all.loc)

#prune based on major allele frequency
#now keep those with major allele frequency > 95% in each pop
calc.max.af<-function(x,y){
	if(max(x[,y])>0)
		m<-max(x[,y]/sum(x[,y]))
	else
		m<-0
	return(m)
}

adult.af<-unlist(lapply(gw.all.loc,calc.max.af,y=3))
names(adult.af)<-names(gw.all.loc)
fem.af<-unlist(lapply(gw.all.loc,calc.max.af,y=4))
names(fem.af)<-names(gw.all.loc)
pop.af<-unlist(lapply(gw.all.loc,calc.max.af,y=5))
names(pop.af)<-names(gw.all.loc)
mal.af<-unlist(lapply(gw.all.loc,calc.max.af,y=6))
names(mal.af)<-names(gw.all.loc)
non.af<-unlist(lapply(gw.all.loc,calc.max.af,y=7))
names(non.af)<-names(gw.all.loc)
juv.af<-unlist(lapply(gw.all.loc,calc.max.af,y=8))
names(juv.af)<-names(gw.all.loc)
prg.af<-unlist(lapply(gw.all.loc,calc.max.af,y=10))
names(prg.af)<-names(gw.all.loc)
mom.af<-unlist(lapply(gw.all.loc,calc.max.af,y=11))
names(mom.af)<-names(gw.all.loc)

adt.af.trim<-adult.af[adult.af<=0.95 & adult.af >= 0.05]
juv.af.trim<-juv.af[juv.af<=0.95 & juv.af >= 0.05]
fem.af.trim<-fem.af[fem.af<=0.95 & fem.af>=0.05]
mal.af.trim<-mal.af[mal.af<=0.95 & mal.af>=0.5]
mom.af.trim<-mom.af[mom.af<=0.95 & mom.af>=0.05]
pop.af.trim<-pop.af[pop.af<=0.95 & pop.af>=0.5]
prg.af.trim<-mal.af[mal.af<=0.95 & mal.af>=0.05]
non.af.trim<-non.af[non.af<=0.95 & non.af>=0.05]
#the comparisons
aj.maf<-names(adt.af.trim[names(adt.af.trim) %in% names(juv.af.trim)])
fm.maf<-names(fem.af.trim[names(fem.af.trim) %in% names(mal.af.trim)])
mo.maf<-names(mom.af.trim[names(mom.af.trim) %in% names(pop.af.trim)])
np.maf<-names(prg.af[names(prg.af.trim) %in% names(non.af.trim)])

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


