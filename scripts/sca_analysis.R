#Author: Sarah P. Flanagan
#Date: 10 August 2015
#Purpose: Analyze SCA data from Stacks

rm(list=ls())
library(gplots)
setwd("E://Docs//SCA")
#############################################################################
#***************************************************************************#
###################################FILES#####################################
#***************************************************************************#
#############################################################################
sca.map<-read.table("E://ubuntushare//SCA//sca_popmap.txt")
sex<-as.character(sca.map[,2])
sex[sex=="FEM"]<-2
sex[sex == "NPM"]<-1
sex[sex == "PRM"]<-1
sex[sex == "OFF"]<-0
sca.ind.info<-data.frame(
	IndID=sub('sample_(\\w{3}\\d{3})_align','\\1',sca.map[,1]), 
	year=rep(2011),locality=rep(1),
	population=sca.map[,2],sex=sex,age=rep(-9),stay=rep(-9),
	recruits=rep(-9),survival=rep(-9))
write.table(sca.ind.info, "E://ubuntushare//SCA//sca.ind.info.txt", 
	col.names=T,row.names=T,eol='\n',sep='\t', quote=F)
popgen.sex.fst.outliers<-read.table("E://Docs//PopGen//sex.outliers.txt", 
	sep='\t',header=T)
popgen.sex.assoc<-as.data.frame(read.table(
	"E://ubuntushare//stacks//populations_sex//sex.plink.assoc",
	header=T))
popgen.sex.map<-read.table(
	"E://ubuntushare//stacks//populations_sex//batch_1.plink.map")

#STACKS SCA
mal.fem.fst<-read.delim(
	"E://ubuntushare//SCA//pops_mal//batch_1.fst_FEM-MAL.tsv")
adt.off.fst<-read.delim(
	"E://ubuntushare//SCA//populations_adt-off//batch_1.fst_ADT-OFF.tsv")
prm.npm.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_NPM-PRM.tsv")
#stacks sca other 
prm.fem.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_FEM-PRM.tsv")
npm.fem.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_FEM-NPM.tsv")
off.fem.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_FEM-OFF.tsv")
prm.off.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_OFF-PRM.tsv")
npm.off.fst<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.fst_NPM-OFF.tsv")
fem.off.fst<-read.delim(
	"E://ubuntushare//SCA//pops_mal//batch_1.fst_FEM-OFF.tsv")
mal.off.fst<-read.delim(
	"E://ubuntushare//SCA//pops_mal//batch_1.fst_MAL-OFF.tsv")

pg.prm.mal<-merge(popgen.sex.fst.outliers, prm.fem.fst, by.x="Chrom",by.y="Chr")
pg.prm.mal<-pg.prm.mal[pg.prm.mal$BP.x == pg.prm.mal$BP.y,] 

#gw_sca output
gw.np.fst<-read.delim(
	"E://GitHub//stacks_processing//gw_sca//gw_sca//gw_sca//sco11_gwsca_Ssco11_within_bool0_bool1.txt")
gw.mf.fst<-read.delim(
	"E://GitHub//stacks_processing//gw_sca//gw_sca//gw_sca//sco11_gwsca_Ssco11_within_sex0_sex1.txt")
gw.oa.fst<-read.delim(
	"E://GitHub//stacks_processing//gw_sca//gw_sca//gw_sca//sco11_gwsca_Ssco11_within_age0_age1.txt")

#############################################################################
#***************************************************************************#
##############################gw_sca FORMATTING##############################
#***************************************************************************#
#############################################################################
#reformat files
sca.ped<-read.delim("E://ubuntushare//SCA//populations//batch_1.plink.ped", 
	skip=1, header=F)
sca.ped$V2<-sub('sample_(\\w{3}\\d{3}.*)_align','\\1',sca.ped$V2)
write.table(sca.ped, quote=F, col.names=F, row.names=F, sep='\t',
	"E://GitHub//stacks_processing//gw_sca//gw_sca//gw_sca//sco11.plink.ped")

sca.ind.info<-read.delim("E://ubuntushare//SCA//sca.ind.info.txt",
	header=T)
#Rename columns
colnames(sca.ind.info)<-c("ID", "year", "locality", "population","sex","age",
	"groups","threshold","bool")
#Recode locality
sca.ind.info$locality[sca.ind.info$population=="PRM"]<-1
sca.ind.info$locality[sca.ind.info$population=="NPM"]<-2
sca.ind.info$locality[sca.ind.info$population=="FEM"]<-3
sca.ind.info$locality[sca.ind.info$population=="OFF"]<-4
#Recode population
sca.ind.info$population<-rep("Ssco11")
#Recode age
sca.ind.info$age[sca.ind.info$locality==4]<-0
sca.ind.info$age[sca.ind.info$locality!=4]<-1
#Recode boolean phenotype
sca.ind.info$bool[sca.ind.info$locality==2]<-0
sca.ind.info$bool[sca.ind.info$locality==1]<-1
#re-write
write.table(sca.ind.info, quote=F, col.names=F, row.names=F, sep='\t',
	"E://GitHub//stacks_processing//gw_sca//gw_sca//gw_sca//sco11.ind.txt")


#############################################################################
#***************************************************************************#
#############################STACKS ANALYSIS#################################
#***************************************************************************#
#############################################################################
#calculate 99% cutoffs
prm.fem.99<-c(mean(prm.fem.fst$Smoothed.Fst)+sd(prm.fem.fst$Smoothed.Fst),
	mean(prm.fem.fst$Smoothed.Fst)-sd(prm.fem.fst$Smoothed.Fst))
npm.fem.99<-c(mean(npm.fem.fst$Smoothed.Fst)+sd(npm.fem.fst$Smoothed.Fst),
	mean(npm.fem.fst$Smoothed.Fst)-sd(npm.fem.fst$Smoothed.Fst))
off.fem.99<-c(mean(off.fem.fst$Smoothed.Fst)+sd(off.fem.fst$Smoothed.Fst),
	mean(off.fem.fst$Smoothed.Fst)-sd(off.fem.fst$Smoothed.Fst))
prm.npm.99<-c(mean(prm.npm.fst$Smoothed.Fst)+sd(prm.npm.fst$Smoothed.Fst),
	mean(prm.npm.fst$Smoothed.Fst)-sd(prm.npm.fst$Smoothed.Fst))
prm.off.99<-c(mean(prm.off.fst$Smoothed.Fst)+sd(prm.off.fst$Smoothed.Fst),
	mean(prm.off.fst$Smoothed.Fst)-sd(prm.off.fst$Smoothed.Fst))
npm.off.99<-c(mean(npm.off.fst$Smoothed.Fst)+sd(npm.off.fst$Smoothed.Fst),
	mean(npm.off.fst$Smoothed.Fst)-sd(npm.off.fst$Smoothed.Fst))
adt.off.99<-c(mean(adt.off.fst$Smoothed.Fst)+sd(adt.off.fst$Smoothed.Fst),
	mean(adt.off.fst$Smoothed.Fst)-sd(adt.off.fst$Smoothed.Fst))
mal.fem.99<-c(mean(mal.fem.fst$Smoothed.Fst)+sd(mal.fem.fst$Smoothed.Fst),
	mean(mal.fem.fst$Smoothed.Fst)-sd(mal.fem.fst$Smoothed.Fst))
fem.off.99<-c(mean(fem.off.fst$Smoothed.Fst)+sd(fem.off.fst$Smoothed.Fst),
	mean(fem.off.fst$Smoothed.Fst)-sd(fem.off.fst$Smoothed.Fst))
mal.off.99<-c(mean(mal.off.fst$Smoothed.Fst)+sd(mal.off.fst$Smoothed.Fst),
	mean(mal.off.fst$Smoothed.Fst)-sd(mal.off.fst$Smoothed.Fst))

prm.fem.sig<-prm.fem.fst[prm.fem.fst$Smoothed.Fst >= prm.fem.99[1] |
	prm.fem.fst$Smoothed.Fst <= prm.fem.99[2],]
npm.fem.sig<-npm.fem.fst[npm.fem.fst$Smoothed.Fst >= npm.fem.99[1] |
	npm.fem.fst$Smoothed.Fst <= npm.fem.99[2],]
off.fem.sig<-off.fem.fst[off.fem.fst$Smoothed.Fst >= off.fem.99[1] |
	off.fem.fst$Smoothed.Fst <= off.fem.99[2],]
prm.npm.sig<-prm.npm.fst[prm.npm.fst$Smoothed.Fst >= prm.npm.99[1] |
	prm.npm.fst$Smoothed.Fst <= prm.npm.99[2],]
npm.off.sig<-npm.off.fst[npm.off.fst$Smoothed.Fst >= npm.off.99[1] |
	npm.off.fst$Smoothed.Fst <= npm.off.99[2],]
adt.off.sig<-adt.off.fst[adt.off.fst$Smoothed.Fst >= adt.off.99[1] |
	adt.off.fst$Smoothed.Fst <= adt.off.99[2],]
mal.fem.sig<-mal.fem.fst[mal.fem.fst$Smoothed.Fst >= mal.fem.99[1] |
	mal.fem.fst$Smoothed.Fst <= mal.fem.99[2],]
fem.off.sig<-fem.off.fst[fem.off.fst$Smoothed.Fst >= fem.off.99[1] |
	fem.off.fst$Smoothed.Fst <= fem.off.99[2],]
mal.off.sig<-mal.off.fst[mal.off.fst$Smoothed.Fst >= mal.off.99[1] |
	mal.off.fst$Smoothed.Fst <= mal.off.99[2],]

#ones that are really diff between males and females
mal.fem<-prm.fem.sig[(prm.fem.sig$Chr %in% npm.fem.sig$Chr &
	prm.fem.sig$BP %in% npm.fem.sig$BP),]
mal.fem<-mal.fem[(mal.fem$Chr %in% mal.fem.sig$Chr &
	mal.fem$BP %in% mal.fem$BP),]
dim(mal.fem[mal.fem$Chr %in% popgen.sex.fst.outliers$Chrom &
	mal.fem$BP %in% popgen.sex.fst.outliers$BP,])#34

#ones that are really diff between adults and offspring
gup.off<-adt.off.sig[adt.off.sig$Chr %in% npm.off.sig$Chr &
	adt.off.sig$BP %in% npm.off.sig$BP,]
gup.off<-gup.off[gup.off$Chr %in% off.fem.sig$Chr &
	gup.off$BP %in% off.fem.sig$BP,]
gup.off<-gup.off[gup.off$Chr %in% fem.off.sig$Chr &
	gup.off$BP %in% fem.off.sig$BP,]
gup.off<-gup.off[gup.off$Chr %in% mal.off.sig$Chr &
	gup.off$BP %in% mal.off.sig$BP,]

##########gw_sca confidence intervals############
gw.np.99<-c(mean(gw.np.fst$WeightedFst)+sd(gw.np.fst$WeightedFst),
	mean(gw.np.fst$WeightedFst)-sd(gw.np.fst$WeightedFst))
gw.mf.99<-c(mean(gw.mf.fst$WeightedFst)+sd(gw.mf.fst$WeightedFst),
	mean(gw.mf.fst$WeightedFst)-sd(gw.mf.fst$WeightedFst))
gw.oa.99<-c(mean(gw.oa.fst$WeightedFst)+sd(gw.oa.fst$WeightedFst),
	mean(gw.oa.fst$WeightedFst)-sd(gw.oa.fst$WeightedFst))

gw.np.out<-gw.np.fst[gw.np.fst$WeightedFst >= gw.np.99[1] |
	gw.np.fst$WeightedFst <= gw.np.99[2],]
gw.np.out$Locus.ID<-sub('(\\d+)_\\d+','\\1',gw.np.out$Locus)
gw.mf.out<-gw.mf.fst[gw.mf.fst$WeightedFst >= gw.mf.99[1] |
	gw.mf.fst$WeightedFst <= gw.mf.99[2],]
gw.mf.out$Locus.ID<-sub('(\\d+)_\\d+','\\1',gw.mf.out$Locus)
gw.oa.out<-gw.oa.fst[gw.oa.fst$WeightedFst >= gw.oa.99[1] |
	gw.oa.fst$WeightedFst <= gw.oa.99[2],]
gw.oa.out$Locus.ID<-sub('(\\d+)_\\d+','\\1',gw.oa.out$Locus)

dim(gw.np.out[gw.np.out$Locus.ID %in% prm.npm.sig$Locus.ID,])
dim(gw.oa.out[gw.oa.out$Locus.ID %in% adt.off.sig$Locus.ID,])
dim(gw.mf.out[gw.mf.out$Locus.ID %in% mal.fem.sig$Locus.ID,])


###########################################################################
#*************************************************************************#
###############################Male offspring##############################
#*************************************************************************#
###########################################################################
prm<-sca.map$V1[sca.map$V2=="PRM"]
prm<-factor(sub('(sample_)(PRM\\d{3}.*)(_align)','\\2',prm))
off<-sca.map$V1[sca.map$V2=="OFF"]
off<-factor(sub('(sample_)(OFF\\d{3}.*)(_align)','\\2',off))

#list of all dads and offspring that were both sequenced. 151 pairs.
dad.off<-c(as.character(prm[sub('\\w{3}(\\d{3}).*','\\1',prm) %in%
	sub('\\w{3}(\\d{3}).*','\\1',off)]),
	as.character(off[sub('\\w{3}(\\d{3}).*','\\1',off) %in%
	sub('\\w{3}(\\d{3}).*','\\1',prm)]))
write.table(dad.off, "dad_off_list.txt", col.names=F, row.names=F,quote=F)

dads<-as.character(prm[sub('\\w{3}(\\d{3}).*','\\1',prm) %in%
	sub('\\w{3}(\\d{3}).*','\\1',off)])
dads[dads=="PRM086-23"]<-"PRM862"
dads<-as.data.frame(cbind(dads,sub('\\w{3}(\\d{3}).*','\\1',dads)))

kids<-as.character(off[sub('\\w{3}(\\d{3}).*','\\1',off) %in%
	sub('\\w{3}(\\d{3}).*','\\1',prm)])
kids[kids=="OFF08623"]<-"OFF862"
kids<-as.data.frame(cbind(kids,sub('\\w{3}(\\d{3}).*','\\1',kids)))

dads.kids<-merge(dads,kids,by.x="V2",by.y="V2")
write.table(dads.kids,"dad.kid.pairs.txt",col.names=F,row.names=F,quote=F)

###########################################################################
#*************************************************************************#
###############################Mom-female Fst##############################
#*************************************************************************#
###########################################################################
mom.fem.sum<-read.delim(
	"E://ubuntushare//SCA//mom_female_fst//fem.mom.fst.summary.txt")
sum.0<-mom.fem.sum[mom.fem.sum$Pop==0,]#females
sum.1<-mom.fem.sum[mom.fem.sum$Pop==1,]#moms
#just use loci that overlap
sum.0<-sum.0[sum.0$Locus %in% sum.1$Locus,]
sum.1<-sum.1[sum.1$Locus %in% sum.0$Locus,]

sum.0.split<-split(sum.0,sum.0$Locus)
sum.1.split<-split(sum.1,sum.1$Locus)

sum.0.hs<-lapply(sum.0.split,function(x){
	hs<-1
	for(i in 1:nrow(x)){
		hs<-1-(x[i,"Freq"]*x[i,"Freq"])
	}
	return(hs)
})
sum.1.hs<-lapply(sum.1.split,function(x){
	hs<-1
	for(i in 1:nrow(x)){
		hs<-1-((x[i,"Count"]/x[i,"PopSize"])*(x[i,"Count"]/x[i,"PopSize"]))
	}
	return(hs)
})
#what if i just split the original thing by locus?
sum.split<-split(mom.fem.sum,mom.fem.sum$Locus)

sum.fst<-lapply(sum.split, function(x){
	if(length(levels(factor(x$Pop))) > 1){
	x.split<-split(x, x$Pop)
	shared.alleles<-factor(x.split[[1]]$Allele[factor(x.split[[1]]$Allele) 
		%in% factor(x.split[[2]]$Allele)])
	private.alleles<-c(as.character(
		factor(x.split[[1]]$Allele[!(factor(x.split[[1]]$Allele) %in% 
			factor(x.split[[2]]$Allele))])), as.character(
		factor(x.split[[2]]$Allele[!(factor(x.split[[2]]$Allele) %in% 
			factor(x.split[[1]]$Allele))])))
	ht<-1
	hs1<-1
	hs2<-1	
	if(length(shared.alleles)>0){
	    for(i in 1:length(shared.alleles)){
		if(length(x.split[[1]][x.split[[1]]$Allele %in%
				shared.alleles[i],"AlleleCount"]) > 0
			& length(x.split[[2]][x.split[[2]]$Allele %in%
				shared.alleles[i],"AlleleCount"]) > 0){
			tot.count<-x.split[[1]][x.split[[1]]$Allele %in% shared.alleles[i],
				"AlleleCount"] +
				x.split[[2]][x.split[[2]]$Allele %in% shared.alleles[i],
				"AlleleCount"]
	    	} 
	    tot.freq<-tot.count/((2*x.split[[1]]$PopSize[1])+
		x.split[[2]]$PopSize[1])
	    ht<-ht-(tot.freq*tot.freq)
	    hs1<-hs1-(x.split[[1]]$Freq[i]*x.split[[2]]$Freq[i])
	    hs2<-hs2-(x.split[[1]]$Freq[i]*x.split[[2]]$Freq[i])
	  }
	} 
	if(length(private.alleles)>0){
	  for(i in 1:length(private.alleles)){
		if(length(x.split[[1]][x.split[[1]]$Allele %in%
			private.alleles[i],"AlleleCount"]) > 0){
			count1<-x.split[[1]][x.split[[1]]$Allele %in%
				private.alleles[i],"AlleleCount"]
			count2<-0
		}
		if(length(x.split[[2]][x.split[[2]]$Allele %in%
				private.alleles[i],"AlleleCount"]) > 0){
			count1<-0
			count2<-x.split[[2]][x.split[[2]]$Allele %in%
				private.alleles[i],"AlleleCount"]
		}
		tot.count<-count1+count2
		tot.freq<-tot.count/((2*x.split[[1]]$PopSize[1])+
			x.split[[2]]$PopSize[1])
	   	ht<-ht-(tot.freq*tot.freq)
	   	hs1<-hs1-((count1/(2*x.split[[1]]$PopSize[1]))*
			(count1/(2*x.split[[1]]$PopSize[1])))
	   	hs2<-hs2-((count2/(x.split[[2]]$PopSize[1]))*
			(count2/(x.split[[2]]$PopSize[1])))
	  }
	}
	hs<-mean(hs1,hs2)
	fst<-(ht-hs)/ht
	return(data.frame(hs1,hs2,hs,ht,fst))}
})
df<-data.frame()
count<-0
for(i in 1:length(sum.fst)){
	if(length(sum.fst[[i]])>0){
		df<-rbind(df,sum.fst[[i]])
		rownames(df)[count]<-names(sum.fst[i])
		count<-count+1
	}
}


mom.fem<-read.delim("E://ubuntushare//SCA//mom_female_fst//fem.mom.fst.txt")
mom.fem.99<-c(mean(mom.fem$Fst)+sd(mom.fem$Fst),
	mean(mom.fem$Fst)-sd(mom.fem$Fst))
mom.fem.out<-mom.fem[mom.fem$Fst >= mom.fem.99[1] |
	mom.fem$Fst <= mom.fem.99[2],]

dim(gw.oa.out[gw.oa.out$Locus.ID %in% mom.fem.out$CatID,])
dim(gw.mf.out[gw.mf.out$Locus.ID %in% mom.fem.out$CatID,])

#create a venn diagram
out.venn<-venn( list("Male-Female"=gw.mf.out$Locus.ID, 
	"Adult-Offspring"=gw.oa.out$Locus.ID, 
	"Pregnant-NonPregnant"=gw.np.out$Locus.ID,
	"Mothers-Females"=mom.fem.out$CatID))

jpeg("venn.comparisons.jpeg", height=7,width=7,units="in", res=300)
par(oma=c(1,1,1,1),cex=0.75)
plot(out.venn, small=0.9)
dev.off()

sca.plink.map<-read.delim(
	"E://ubuntushare//SCA//populations//batch_1.plink.map", skip=1, header=F)
sca.plink.map$loc<-sub('(\\d+)_\\d+','\\1',sca.plink.map$V2)
colnames(sca.plink.map)<-c("chr","snp","dist","bp","loc")
gw.np.fst<-merge(gw.np.fst, sca.plink.map,by.x="Locus", by.y="snp")
gw.oa.fst<-merge(gw.oa.fst, sca.plink.map,by.x="Locus",by.y="snp")
gw.mf.fst<-merge(gw.mf.fst,sca.plink.map,by.x="Locus",by.y="snp")
mom.fem.fst<-merge(mom.fem, sca.plink.map, by.x="CatID",by.y="loc")
#used haplotypes, not snps, so need to remove duplicates & keep one bp per locus
mom.fem.fst<-mom.fem.fst[!(duplicated(mom.fem.fst$CatID)),]
#there are 10060 catids that aren't in the sca.plink.map...??


jpeg("fst.outliers.jpeg", height=10,width=7.5,units="in",res=300)
par(mfrow=c(4,1),oma=c(1,1,0,0),mar=c(0,1,1,0),cex=2,mgp=c(3,0.5,0))
plot.fsts(gw.np.fst,gw.np.99,fst.name="WeightedFst",bp.name="bp")
plot.fsts(gw.mf.fst,gw.mf.99,fst.name="WeightedFst",bp.name="bp")
plot.fsts(gw.oa.fst,gw.oa.99,fst.name="WeightedFst",bp.name="bp")
#plot.fsts(mom.fem.fst, mom.fem.99,fst.name="Fst",bp.name="bp",chrom.name="chr")
chrom.name="chr"
bp.name="bp"
fst.name="Fst"
all.scaff<-split(mom.fem.fst, mom.fem.fst[,"chr"])
map.split<-split(sca.plink.map, sca.plink.map$chr)
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	for(i in 1:length(map.split)){
		new.max<-last.max+round(max(map.split[[i]][,bp.name]), -2)
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	#change BP to plot
	for(i in 1:length(map.split)){
		map.split[[i]][,bp.name]<-
			map.split[[i]][,bp.name]+addition.values[i]
	}
	x.max<-max(addition.values)
	x.min<-min(map.split[[1]][,bp.name])
	y.max<-max(mom.fem.fst[,fst.name])+0.1*max(mom.fem.fst[,fst.name])
	y.min<-min(mom.fem.fst[,fst.name])-0.1*min(mom.fem.fst[,fst.name])
	if(min(mom.fem.fst[,fst.name]) < 0) {
		y.min<-min(mom.fem.fst[,fst.name]) - 0.1*min(mom.fem.fst[,fst.name])
	} else {
		y.min<-0
	}

	plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), 
		ylim=c(y.min, y.max), 
		bty="n",type="n",	axes=F, xlab="", ylab="")
	for(i in 1:nrow(rect.xs)){
		if(i%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
			col=rect.color, border=NA)
	}
	for(i in 1:length(all.scaff)){
		map.bp<-map.split[[i]][map.split[[i]]$loc %in% 
			all.scaff[[i]]$CatID,]
		map.bp<-map.bp[!duplicated(map.bp$loc),]
		plot.genome.wide(map.bp$bp, 
			all.scaff[[i]][,fst.name],plot.rect=FALSE,
			y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col="grey53",
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
		temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] >= ci.dat[1],]
		temp.bp<-map.split[[i]][map.split[[i]]$loc %in% 
			temp.sig$CatID,]
		temp.bp<-temp.bp[!duplicated(temp.bp$loc),]
		points(temp.bp$bp, temp.sig[,fst.name], 
			col="red", pch=19, cex=0.5)
		temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] <= ci.dat[2],]
		temp.bp<-map.split[[i]][map.split[[i]]$loc %in% 
			temp.sig$CatID,]
				temp.bp<-temp.bp[!duplicated(temp.bp$loc),]
		points(temp.bp$bp, temp.sig[,fst.name], 
			col="yellow", pch=19, cex=0.5)
	}
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		labels=seq(round(y.min,2),round(y.max,2),
			round((y.max-y.min)/2, digits=2)),
		las=1,tck = -0.05, xlab="", ylab="", cex.axis=0.5)
	mtext(side=1, "Genomic Location", outer = FALSE, cex=1,line=-0.5)
	mtext(side=2, "Global Fst", outer=FALSE, line=1.2,cex=1)
dev.off()



###########################################################################
#*************************************************************************#
###############################gwSCA haplotypes############################
#*************************************************************************#
###########################################################################
setwd("E:\\ubuntushare\\SCA\\gwsca_haplotypes")
gw.fst<-read.delim("gwsca_fsts.txt")
gw.sum<-read.delim("gwsca_summary.txt")
gw.fst.lump<-read.delim("gwsca_fsts_lumped.txt")
gw.alleles<-read.delim("gwsca_haplotypes_alleles.txt")
gw.cat<-read.delim("E:\\ubuntushare\\SCA\\batch_1.catalog.tags.tsv",header=F)

#Summarize number of alleles
summary(gw.sum[(gw.sum$LocusName %in% aj.loci) & 
	(gw.sum$Pop=="ADULT" | gw.sum$Pop=="JUVIE")
	,"NumAlleles"])
summary(gw.sum[(gw.sum$LocusName %in% fm.loci) &
	(gw.sum$Pop=="MAL" | gw.sum$Pop=="FEM")
	,"NumAlleles"])
summary(gw.sum[(gw.sum$LocusName %in% mo.loci) &
	(gw.sum$Pop=="MOM" | gw.sum$Pop=="POP")
	,"NumAlleles"])
summary(gw.sum[(gw.sum$LocusName %in% np.loci) &
	(gw.sum$Pop=="PREGGER" | gw.sum$Pop=="NONPREG")
	,"NumAlleles"])

###########PRUNING THE DATA##########
#create different structures
gw.all.loc<-split(gw.alleles, gw.alleles$Locus)
gw.all.sum<-lapply(gw.all.loc, function(x) colSums(x[,3:11]))
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
prg.af<-unlist(lapply(gw.all.loc,calc.max.af,y=10))
names(prg.af)<-names(gw.all.loc)
mom.af<-unlist(lapply(gw.all.loc,calc.max.af,y=11))
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
write.table(all.pruned, "pruned_alleles_gwsca.txt",header=T,row.names=F,
	quote=F,sep='\t')
###############PLOTTING################
#generate plotting structure
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


#Calculate Ht
gw.all.loc.aj<-gw.all.loc[names(gw.all.loc) %in% aj.prune]
calc.freq<-function(df, index){
	df.new<-df[,index]/sum(df[,index])
	return(df.new)
}
gw.freq.aj<-lapply(gw.all.loc.aj, function(x){
	df.new<-cbind(x[,c(1:2)],calc.freq(x,3),calc.freq(x,8))
})
calc.ht<-function(df,grp1,grp2){
	ht<-1
	for(i in 1:nrow(df)){
		avg<-(df[i,grp1]+df[i,grp2])/2
		ht<-ht-(avg*avg)
	}
	return(ht)
}
gw.aj.ht<-lapply(gw.freq.aj, calc.ht, grp1=3,grp2=4)
aj.ht<-data.frame(cbind(names(gw.all.loc.aj),as.numeric(unlist(gw.aj.ht))),
	stringsAsFactors=F)
colnames(aj.ht)<-c("Locus","Ht")
aj.ht$Ht<-as.numeric(aj.ht$Ht)

#is Ht evenly distributed throughout the genome? YES
aj.sum<-merge(aj.ht, aj.plot, by.x="Locus",by.y="Locus")
aj.sum$Ht<-as.numeric(aj.sum$Ht)
ci.dat<-c(mean(aj.sum$Ht)+(2.57583*sd(aj.sum$Ht)),
	mean(aj.sum$Ht)-(2.57583*sd(aj.sum$Ht)))
plot.fsts(aj.sum, ci.dat,fst.name="Ht", chrom.name="Chrom")
legend("top","Adult-Juveinle Ht", cex=0.5,bty="n")

#Calculate Fst for grouped minor alleles
calc.lumped.fst<-function(af.1, af.2, loci){
	#only use loci in both sets
	af1<-af.1[names(af.1) %in% names(af.2)]
	af2<-af.2[names(af.2) %in% names(af.1)]
	fst.df<-data.frame(p1=af1[names(af1) %in% loci], 
		q1=1-af1[names(af1) %in% loci],
		p2=af2[names(af2) %in% loci], 
		q2=1-af2[names(af2) %in% loci])
	fst.df$locus<-names(af1[names(af1) %in% loci])
	fst.df$hs1<-2*fst.df$p1*fst.df$q1
	fst.df$hs2<-2*fst.df$p2*fst.df$q2
	fst.df$hs<-(fst.df$hs1+fst.df$hs2)/2
	fst.df$ht<-2*((fst.df$p1+fst.df$p2)/2)*
		((fst.df$q1+fst.df$q2)/2)
	fst.df$fst<-(fst.df$ht-fst.df$hs)/fst.df$ht
	return(fst.df)
}

#af1=adult.af, loci=aj.plot$Locus, af2=juv.af, fst.df=aj.fst.grp
aj.fst.grp<-calc.lumped.fst(adt.prune, juv.prune, aj.plot$Locus)
aj.fst.grp<-merge(aj.fst.grp, aj.plot, by.x="locus", by.y="Locus")
fm.fst.grp<-calc.lumped.fst(fem.prune, mal.prune, fm.plot$Locus)
fm.fst.grp<-merge(fm.fst.grp, fm.plot, by.x="locus", by.y="Locus")
mo.fst.grp<-calc.lumped.fst(mom.prune, pop.prune, mo.plot$Locus)
mo.fst.grp<-merge(mo.fst.grp, mo.plot, by.x="locus", by.y="Locus")
np.fst.grp<-calc.lumped.fst(prg.prune, non.prune, np.plot$Locus)
np.fst.grp<-merge(np.fst.grp, np.plot, by.x="locus", by.y="Locus")

aj.grp.ci<-c(mean(aj.fst.grp$fst)+2.57583*sd(aj.fst.grp$fst),
	mean(aj.fst.grp$fst)-2.57583*sd(aj.fst.grp$fst))
fm.grp.ci<-c(mean(fm.fst.grp$fst)+2.57583*sd(fm.fst.grp$fst),
	mean(fm.fst.grp$fst)-(2.57583*sd(fm.fst.grp$fst)))
mo.grp.ci<-c(mean(mo.fst.grp$fst)+(2.57583*sd(mo.fst.grp$fst)),
	mean(mo.fst.grp$fst)-(2.57583*sd(mo.fst.grp$fst)))
np.grp.ci<-c(mean(np.fst.grp$fst)+(2.57583*sd(np.fst.grp$fst)),
	mean(np.fst.grp$fst)-(2.57583*sd(np.fst.grp$fst)))

png("fst.haps.lumped.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(4,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0),cex=1.5)
#all used to be gw.plot
aj.fst.grp$Chrom<-factor(aj.fst.grp$Chrom)
plot.fsts(aj.fst.grp, ci.dat=aj.grp.ci,fst.name="fst", chrom.name="Chrom")
abline(h=aj.null99[1],lwd=2,col="purple")
legend("top","Adult-Juvenile", cex=0.5,bty="n")
fm.fst.grp$Chrom<-factor(fm.fst.grp$Chrom)
plot.fsts(fm.fst.grp, ci.dat=fm.grp.ci,fst.name="fst", chrom.name="Chrom")
abline(h=fm.null99[1],lwd=2,col="purple")
legend("top","Male-Female", cex=0.5,bty="n")
mo.fst.grp$Chrom<-factor(mo.fst.grp$Chrom)
plot.fsts(mo.fst.grp, ci.dat=mo.grp.ci,fst.name="fst", chrom.name="Chrom")
abline(h=mo.null99[1],lwd=2,col="purple")
legend("top","Mothers-Females", cex=0.5,bty="n")
np.fst.grp$Chrom<-factor(np.fst.grp$Chrom)
plot.fsts(np.fst.grp, ci.dat=np.grp.ci,fst.name="fst", chrom.name="Chrom")
abline(h=np.null99[1],lwd=2,col="purple")
legend("top","Pregnant-NonPregnant", cex=0.5,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

#*************************************************************************#
###########################gwSCA simulation model##########################
#*************************************************************************#
aj.plot<-read.table("E://ubuntushare//SCA//gwsca_haplotypes//aj.plot")
fm.plot<-read.table("E://ubuntushare//SCA//gwsca_haplotypes//fm.plot")
gw.plot<-read.table("E://ubuntushare//SCA//gwsca_haplotypes//gw.plot")

#rand AFS sim
rand.sim<-read.delim("E://ubuntushare//SCA//sca_simulation//randdist.ss0.test.fst_out.txt")
rand.sim$Chrom<-c(rep("0",max(rand.sim$Locus)+1),
	rep("1",max(rand.sim$Locus)+1),rep("2",max(rand.sim$Locus)+1),
	rep("3",max(rand.sim$Locus)+1))
rand.ao<-read.delim("sca_simulation//randdist.ss0.testao.txt")

aj.null99<-c(mean(rand.sim$AOFst)+2.57583*sd(rand.sim$AOFst),
	mean(rand.sim$AOFst)-2.57583*sd(rand.sim$AOFst))
fm.null99<-c(mean(rand.sim$MFFst)+2.57583*sd(rand.sim$MFFst),
	mean(rand.sim$MFFst)-2.57583*sd(rand.sim$MFFst))
mo.null99<-c(mean(rand.sim$MDFst)+2.57583*sd(rand.sim$MDFst),
	mean(rand.sim$MDFst)-2.57583*sd(rand.sim$MDFst))

plot.fsts(rand.sim, ci.dat=aj.null99,fst.name="AOFst",chrom.name="Chrom",bp.name="Locus")
plot.fsts(aj.plot, ci.dat=aj.null99,fst.name="Adult.Juvie", chrom.name="Chrom")

#known AFS sim
known.sim<-read.delim("E://ubuntushare//SCA//sca_simulation//knowndist.ss0.test.fst_out.txt")
known.sim$Chrom<-c(rep("0",max(known.sim$Locus)+1),
	rep("1",max(known.sim$Locus)+1),rep("2",max(known.sim$Locus)+1),
	rep("3",max(known.sim$Locus)+1))
aj.null99<-c(mean(known.sim$AOFst)+2.57583*sd(known.sim$AOFst),
	mean(known.sim$AOFst)-2.57583*sd(known.sim$AOFst))
fm.null99<-c(mean(known.sim$MFFst)+2.57583*sd(known.sim$MFFst),
	mean(known.sim$MFFst)-2.57583*sd(known.sim$MFFst))
mo.null99<-c(mean(known.sim$MDFst)+2.57583*sd(known.sim$MDFst),
	mean(known.sim$MDFst)-2.57583*sd(known.sim$MDFst))

#prune it down to ones in correct AFS
known.sim$AvgAF<-(known.sim$AdultAF+known.sim$OffAF+known.sim$MaleAF
	+ known.sim$FemAF + known.sim$DadAF)/5
sim.keep<-known.sim[known.sim$AvgAF > 0.05 & known.sim$AvgAF < 0.95,]
aj.null99<-c(mean(sim.keep$AOFst)+2.57583*sd(sim.keep$AOFst),
	mean(sim.keep$AOFst)-2.57583*sd(sim.keep$AOFst))
fm.null99<-c(mean(sim.keep$MFFst)+2.57583*sd(sim.keep$MFFst),
	mean(sim.keep$MFFst)-2.57583*sd(sim.keep$MFFst))
mo.null99<-c(mean(sim.keep$MDFst)+2.57583*sd(sim.keep$MDFst),
	mean(sim.keep$MDFst)-2.57583*sd(sim.keep$MDFst))

plot.fsts(known.sim, ci.dat=aj.null99,fst.name="AOFst",chrom.name="Chrom",bp.name="Locus")
plot.fsts(aj.fst.grp, ci.dat=aj.null99,fst.name="Adult.Juvie",chrom.name="Chrom",bp.name="BP")


png("E://Docs//SCA//fst.haps.knownafsmodel.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0),cex=1.5)
plot.fsts(aj.plot, ci.dat=aj.null99,fst.name="Adult.Juvie", chrom.name="Chrom")
legend("top","Adult-Juvenile", cex=1,bty="n")
plot.fsts(fm.plot, ci.dat=fm.null99,fst.name="Fem.Mal", chrom.name="Chrom")
legend("top","Males-Females", cex=1,bty="n")
plot.fsts(mo.plot, ci.dat=mo.null99,fst.name="Fem.Mom", chrom.name="Chrom")
legend("top","Mothers-Females", cex=1,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()


##testing
#the distribution
distribution<-read.table("E://ubuntushare//SCA//sca_simulation//test_empirical_randnums//test_distribution.txt")
hist(distribution$V1)
#the model
known.sim<-read.delim("E://ubuntushare//SCA//sca_simulation//knowndist.ss0.test.fst_out.txt")
known.sim$avg.maf<-(known.sim$AdultAF+known.sim$OffAF+known.sim$MaleAF+ known.sim$FemAF + known.sim$DadAF)/5
hist(known.sim$avg.maf)
known.sim$Chrom<-c(rep("0",max(known.sim$Locus)+1),
	rep("1",max(known.sim$Locus)+1),rep("2",max(known.sim$Locus)+1),
	rep("3",max(known.sim$Locus)+1))
aj.null99<-c(mean(known.sim$AOFst)+2.57583*sd(known.sim$AOFst),
	mean(known.sim$AOFst)-2.57583*sd(known.sim$AOFst))
fm.null99<-c(mean(known.sim$MFFst)+2.57583*sd(known.sim$MFFst),
	mean(known.sim$MFFst)-2.57583*sd(known.sim$MFFst))
mo.null99<-c(mean(known.sim$MDFst)+2.57583*sd(known.sim$MDFst),
	mean(known.sim$MDFst)-2.57583*sd(known.sim$MDFst))

aj.known<-known.sim[known.sim$AOFst >0,]
aj.null99<-c(mean(aj.known$AOFst)+2.57583*sd(aj.known$AOFst),
	mean(aj.known$AOFst)-2.57583*sd(aj.known$AOFst))
fm.known<-known.sim[known.sim$MFFst>0,]
fm.null99<-c(mean(fm.known$MFFst)+2.57583*sd(fm.known$MFFst),
	mean(fm.known$MFFst)-2.57583*sd(fm.known$MFFst))
mo.known<-known.sim[known.sim$MDFst>0,]
mo.null99<-c(mean(mo.known$MDFst)+2.57583*sd(mo.known$MDFst),
	mean(mo.known$MDFst)-2.57583*sd(mo.known$MDFst))


plot.fsts(known.sim, ci.dat=aj.null99,fst.name="AOFst",chrom.name="Chrom",bp.name="Locus")
plot.fsts(aj.plot, ci.dat=aj.null99,fst.name="Adult.Juvie", chrom.name="Chrom")
plot.fsts(fm.plot, ci.dat=fm.null99,fst.name="Fem.Mal", chrom.name="Chrom")
plot.fsts(mo.plot, ci.dat=mo.null99,fst.name="Fem.Mom", chrom.name="Chrom")



#*************************************************************************#
###############################gwSCA null model############################
#*************************************************************************#
setwd("E:\\ubuntushare\\SCA\\null_model_numerical")
aj.null<-read.delim("a440_b160_sampledpops.txt")
fm.null<-read.delim("a244_b196_sampledpops.txt")
mo.null<-read.delim("a157_b87_sampledpops.txt")
np.null<-read.delim("a188_b8_sampledpops.txt")

aj.null99<-c(mean(aj.null$WrightsFst)+2.57583*sd(aj.null$WrightsFst),
	mean(aj.null$WrightsFst)-2.57583*sd(aj.null$WrightsFst))
fm.null99<-c(mean(fm.null$WrightsFst)+2.57583*sd(fm.null$WrightsFst),
	mean(fm.null$WrightsFst)-2.57583*sd(fm.null$WrightsFst))
mo.null99<-c(mean(mo.null$WrightsFst)+2.57583*sd(mo.null$WrightsFst),
	mean(mo.null$WrightsFst)-2.57583*sd(mo.null$WrightsFst))
np.null99<-c(mean(np.null$WrightsFst)+2.57583*sd(np.null$WrightsFst),
	mean(np.null$WrightsFst)-2.57583*sd(np.null$WrightsFst))



png("fst.haps.null.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(4,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0))
#all used to be gw.plot
aj.plot$Chrom<-factor(aj.plot$Chrom)
plot.fsts(aj.plot, ci.dat=aj.null99,fst.name="Adult.Juvie", chrom.name="Chrom")
legend("top","Adult-Juvenile", cex=0.5,bty="n")
fm.plot$Chrom<-factor(fm.plot$Chrom)
plot.fsts(fm.plot, ci.dat=fm.null99,fst.name="Fem.Mal", chrom.name="Chrom")
legend("top","Male-Female", cex=0.5,bty="n")
mo.plot$Chrom<-factor(mo.plot$Chrom)
plot.fsts(mo.plot, ci.dat=mo.null99,fst.name="Fem.Mom", chrom.name="Chrom")
legend("top","Mothers-Females", cex=0.5,bty="n")
np.plot$Chrom<-factor(np.plot$Chrom)
plot.fsts(np.plot, ci.dat=np.null99,fst.name="Nonpreg.Pregger", chrom.name="Chrom")
legend("top","Pregnant-NonPregnat", cex=0.5,bty="n")
dev.off()




aj.null<-read.delim("a440_b160_d1_sampledpops.txt")

ibs.null<-read.delim("E://ubuntushare//SCA//null_model_numerical//sca_sim_output//null_10reps_440a_160j_6alleles_plot_Fsts.txt")
ibs.null<-read.delim("E:\\ubuntushare\\simulation_model\\simulation_model\\simulation_model\\null_1rep_440a_160j_6alleles_plot_Fsts.txt")
ibs.sum<-read.delim("E:\\ubuntushare\\simulation_model\\simulation_model\\simulation_model\\null_1rep_440a_160j_6alleles_SummStats.txt")
ibs.fst<-read.table("E:\\ubuntushare\\simulation_model\\simulation_model\\simulation_model\\null_1rep_440a_160j_6alleles_Fsts.txt",skip=4,header=T)
ibs.null99<-mean(ibs.null$Smooth99)
plot.fsts(aj.plot, ci.dat=c(ibs.null99,0),fst.name="Adult.Juvie", chrom.name="Chrom")
plot.fsts(aj.fst.grp, ci.dat=c(ibs.null99,0),fst.name="fst", chrom.name="Chrom")



png("E:\\Docs\\SCA\\aj.fst.stacks.null.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(2,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0))
plot.fsts(adt.off.fst, ci.dat=c(aj.null99[1],-5),fst.name="Smoothed.Fst", chrom.name="Chr")
legend("top","Numerical Analysis", cex=0.75,bty="n")
plot.fsts(adt.off.fst, ci.dat=ibs.null99,fst.name="Smoothed.Fst", chrom.name="Chr")
legend("top","Individual-Based Simulations", cex=0.75,bty="n")
dev.off()

png("E:\\Docs\\SCA\\mf.fst.null.png",height=150,width=300,units="mm",res=300)
plot.fsts(gw.plot, ci.dat=fm.null99,fst.name="Fem.Mom", chrom.name="Chrom")
dev.off()

plot.fsts(gw.plot, ci.dat=c(10,-5),fst.name="Fem.Mom", chrom.name="Chrom")


######from Dec 21-23
##figuring out coverage for monnahan sca
vcf<-read.table(
	"E://ubuntushare//SCA//gatk//filter_vcf//genotype_output_filtered.vcf",
	skip=3603, header=T)

vcf.summ<-read.table(
	"E://ubuntushare//SCA//gatk//filter_vcf//filtered_summary.txt", 
	header=T)

stacks.coverage<-read.table(
	"E://empty_vs_project//coverage",
	header=T)

locus.info<-data.frame(locus=gw.cat$V3, scaffold=gw.cat$V4, bp=gw.cat$V5)

stacks.info<-merge(stacks.coverage, locus.info, by.x="ID", by.y="locus")

info<-merge(stacks.info, vcf.summ, by.x=c("scaffold","bp"),by.y=c("Scaffold","Pos"))




###########################################################################
#*************************************************************************#
################################Merge VCFs#################################
#*************************************************************************#
###########################################################################
mal.vcf<-read.table("E://ubuntushare//SCA//batch_1.mal.vcf")
fem.vcf<-read.table("E://ubuntushare//SCA//batch_1.fem.vcf")
off.vcf<-read.table("E://ubuntushare//SCA//batch_1.off.vcf")

mal.fem.vcf<-merge(mal.vcf,fem.vcf,c("V1","V2","V3"))
merge.vcf<-merge(mal.fem.vcf,off.vcf,c("V1","V2","V3"))
vcf.loc<-merge.vcf[,1:3]
colnames(vcf.loc)<-c("Chr","Pos","ID")

mal.vcf<-merge(mal.vcf,vcf.loc,c("V1","V2","V3"))
fem.vcf<-merge(fem.vcf,vcf.loc,c("V1","V2","V3"))
off.vcf<-merge(off.vcf,vcf.loc,c("V1","V2","V3"))

write.table(mal.vcf,"E://ubuntushare//SCA//mal.vcf", row.names=F,col.names=F,
	quote=F, sep='\t')
write.table(fem.vcf,"E://ubuntushare//SCA//fem.vcf", row.names=F,col.names=F,
	quote=F, sep='\t')
write.table(off.vcf,"E://ubuntushare//SCA//off.vcf", row.names=F,col.names=F,
	quote=F, sep='\t')

#create a whitelist for populations
write.table(vcf.loc$ID,"E://ubuntushare//SCA//shared_loci.txt",row.names=F,
	col.names=F, quote=F,sep = '\t', eol='\n')

sca.tag<-read.table("E://ubuntushare//SCA//batch_1.catalog.tags.tsv")
sca.snp<-read.table("E://ubuntushare//SCA//batch_1.catalog.snps.tsv")
sca.ref<-merge(sca.tag[,3:5], sca.snp[,3:4],by="V3")
colnames(sca.ref)<-c("LocusID","Chr","BP","Col")
sca.vcf.white<-merge(sca.ref,vcf.loc,by.x=c("Chr","BP"),by.y=c("Chr","Pos"))


###########################################################################
#*************************************************************************#
#############################MONNAHAN ANALYSIS#############################
#*************************************************************************#
###########################################################################
#plotting het_v_depth output
hd<-read.table("E://ubuntushare//SCA//monnahan//het_v_depth.txt", sep="\t", header=T)

hd$prop.het<-hd$called.het/(hd$called.het+hd$called.homo)
hd$prop.hom<-hd$called.homo/(hd$called.het+hd$called.homo)
plot(hd$read.depth, hd$prop.het)
points(hd$read.depth, hd$prop.hom, col="blue")

#compare loci from gatk to those in sca
setwd("E://ubuntushare//SCA//results//")

sumstats<-read.table("stacks/batch_1.sumstats.tsv")
inform.snps<-sumstats[,2:5]
colnames(inform.snps)<-c("LocID","Chr","BP","Col")
inform.snps$SNPinfo<-paste(inform.snps$Chr,"_",inform.snps$BP,sep="")
rm(sumstats)

gvcf<-read.table("monnahan/out.vcf")
gvcf.info<-data.frame(Chr=gvcf$V1, BP=gvcf$V2)
gvcf.info$SNPinfo<-paste(gvcf.info$Chr,"_",gvcf.info$BP,sep="")

shared<-merge(gvcf.info, inform.snps, by="SNPinfo")
write.table(shared,"stacks.gatk.shared.txt",quote=F,
	row.names=F,col.names=F,sep="\t")


###########################################################################
#*************************************************************************#
#########################BIALLELIC GWSCA WITH VCF##########################
#*************************************************************************#
###########################################################################
vcf<-read.delim("E://ubuntushare//SCA//monnahan//batch_1.vcf", skip=9)
dad.kid.pairs<-read.table("E://Docs//SCA//dad.kid.pairs.txt", header=F,
	stringsAsFactors=F)
dad.kid.pairs<-dad.kid.pairs[,-1]
colnames(dad.kid.pairs)<-c("Dads","Kids")

colnames(vcf)[colnames(vcf)=="sample_PRM086.23_align"]<-"sample_PRM8623_align"
colnames(vcf)[colnames(vcf)=="sample_OFF08623_align"]<-"sample_OFF8623_align"
colnames(vcf)[colnames(vcf)=="sample_PRM035.1_align"]<-"sample_PRM035_align"
dad.kid.pairs$Dads[dad.kid.pairs$Dads=="PRM035-1"]<-"PRM035"

dads<-cbind(vcf[,1:9],vcf[,
	sub('sample_(\\w{3}\\d{3})*_align','\\1',colnames(vcf))
	 %in% dad.kid.pairs$Dads])
kids<-cbind(vcf[,1:9],vcf[,
	sub('sample_(\\w{3}\\d{3}).*_align','\\1',colnames(vcf))
	 %in% dad.kid.pairs$Kids])

keep.columns<-c(colnames(dads)[10:ncol(dads)],colnames(kids)[10:ncol(kids)])
write.csv(keep.columns,"E://ubuntushare//SCA//keep.columns.csv",
	row.names=F)


























