
###########################################################################
#*************************************************************************#
###############################gwSCA haplotypes############################
#*************************************************************************#
###########################################################################
source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/haplotypes")
hap.fst<-read.delim("gwsca_hap_fsts.txt")
hap.sum<-read.delim("gwsca_hap_summary.txt")
gw.fst.lump<-read.delim("gwsca_fsts_lumped.txt")
gw.alleles<-read.delim("gwsca_haplotypes_alleles.txt")
gw.cat<-read.delim("../stacks/batch_1.catalog.tags.tsv",header=F)

#Summarize number of alleles
summary(hap.sum[(hap.sum$LocusName %in% aj.loci) & 
	(hap.sum$Pop=="ADULT" | hap.sum$Pop=="JUVIE")
	,"NumAlleles"])
summary(hap.sum[(hap.sum$LocusName %in% fm.loci) &
	(hap.sum$Pop=="MAL" | hap.sum$Pop=="FEM")
	,"NumAlleles"])
summary(hap.sum[(hap.sum$LocusName %in% mo.loci) &
	(hap.sum$Pop=="MOM" | hap.sum$Pop=="POP")
	,"NumAlleles"])
summary(hap.sum[(hap.sum$LocusName %in% np.loci) &
	(hap.sum$Pop=="PREGGER" | hap.sum$Pop=="NONPREG")
	,"NumAlleles"])

###########PRUNING THE DATA##########
#create different structures
gw.all.loc<-split(gw.alleles, gw.alleles$Locus)
gw.all.sum<-lapply(gw.all.loc, function(x) colSums(x[,3:10]))
gw.loc.count<-data.frame(matrix(unlist(gw.all.sum), nrow=length(gw.all.sum), byrow=T),stringsAsFactors=FALSE)
colnames(gw.loc.count)<-names(gw.all.sum[[1]])
rownames(gw.loc.count)<-names(gw.all.loc)

#prune based on major allele frequency
#now keep those with major allele frequency > 95% in each pop
#based on: max(gw.all.loc[[1]]$ADULTCount/sum(gw.all.loc[[1]]$ADULTCount))
calc.max.af<-function(x,y){
	if(max(x[,y])>0)
		m<-max(x[,y]/sum(x[,y]))
	else
		m<-0
	return(m)
}
#calculate maximum allele frequencies
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
prg.af<-unlist(lapply(gw.all.loc,calc.max.af,y=9))
names(prg.af)<-names(gw.all.loc)
mom.af<-unlist(lapply(gw.all.loc,calc.max.af,y=10))
names(mom.af)<-names(gw.all.loc)

#prune allele frequencies (this method works!)
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
mo.maf<-names(mom.af.trim[names(mom.af.trim) %in% names(fem.af.trim)])
np.maf<-names(prg.af[names(prg.af.trim) %in% names(non.af.trim)])

##NOW prune based on representation in the groups
adt.grp.trim<-rownames(gw.loc.count[gw.loc.count$ADULTCount>216,])
juv.grp.trim<-rownames(gw.loc.count[gw.loc.count$JUVIECount>159,])
fem.grp.trim<-rownames(gw.loc.count[gw.loc.count$FEMCount>57,])
mal.grp.trim<-rownames(gw.loc.count[gw.loc.count$MALCount>159,])
mom.grp.trim<-rownames(gw.loc.count[gw.loc.count$MOMCount>130,])
pop.grp.trim<-rownames(gw.loc.count[gw.loc.count$POPCount>57,])
prg.grp.trim<-rownames(gw.loc.count[gw.loc.count$PREGGERCount>159,])
non.grp.trim<-rownames(gw.loc.count[gw.loc.count$NONPREGCount>8,])
#the comparisons
aj.loci<-rownames(gw.loc.count[gw.loc.count$ADULTCount>216 & 
	gw.loc.count$JUVIECount>159,])
fm.loci<-rownames(gw.loc.count[gw.loc.count$FEMCount>57 & 
	gw.loc.count$MALCount>159,])
mo.loci<-rownames(gw.loc.count[gw.loc.count$FEMCount>57 & 
	gw.loc.count$MOMCount>25,])
np.loci<-rownames(gw.loc.count[gw.loc.count$PREGGERCount>159 & 
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
write.table(all.pruned, "pruned_alleles_gwsca.txt",header=T,row.names=F,
	quote=F,sep='\t')
###############PLOTTING################
#generate plotting structure
gw.plot<-data.frame(Locus=hap.fst$Locus,Adult.Juvie=hap.fst$ADULT.JUVIE, 
	Fem.Mal=hap.fst$FEM.MAL, Fem.Mom=hap.fst$FEM.MOM, 
	Nonpreg.Pregger=hap.fst$NONPREG.PREGGER)
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

aj.plot$Chrom<-factor(aj.plot$Chrom)
fm.plot$Chrom<-factor(fm.plot$Chrom)
mo.plot$Chrom<-factor(mo.plot$Chrom)
np.plot$Chrom<-factor(np.plot$Chrom)
#plot
png("fst.hap.filtered.distci.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(4,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(aj.plot, ci.dat=aj.ci,fst.name="Adult.Juvie", chrom.name="Chrom"
	, axis.size=0.75)
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.plot, ci.dat=fm.ci,fst.name="Fem.Mal", chrom.name="Chrom"
	, axis.size=0.75)
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(mo.plot, ci.dat=mo.ci,fst.name="Fem.Mom", chrom.name="Chrom"
	, axis.size=0.75)
legend("top","Mothers-Females", cex=0.75,bty="n")
plot.fsts(np.plot, ci.dat=np.ci,fst.name="Nonpreg.Pregger",chrom.name="Chrom",
	axis.size=0.75)
legend("top","Pregnant-NonPregnant", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

