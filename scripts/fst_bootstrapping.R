setwd("C:/Users/Sarah/OneDrive")

fst<-read.delim("gwsca_fsts.txt")
fst$Locus<-paste(fst$Chrom, fst$Pos,sep=".")
sum<-read.delim("gwsca_summary.txt")
sum$Locus<-paste(sum$Chrom,sum$Pos,sep=".")

fst<-fst[,c("Chrom","Pos","Locus","ADULT.JUVIE","FEM.MAL","FEM.MOM",
	"NONPREG.PREGGER","MAL.JUVIE")]
fst.split<-split(fst,fst$Chrom)

######PRUNE#######
#remove any that are not polymorphic
sum<-sum[!is.na(sum$AA),]
#test for HWE
sum$AAexp<-sum$Allele1Freq*sum$Allele1Freq
sum$aaexp<-sum$Allele2Freq*sum$Allele2Freq
sum$Aaexp<-1-sum$aaexp-sum$AAexp
sum$chi<-(((sum$AA-sum$AAexp)^2)/sum$AAexp)+
	(((sum$Aa-sum$Aaexp)^2)/sum$Aaexp)+
	(((sum$aa-sum$aaexp)^2)/sum$aaexp)
sum$chi.result<-1-pchisq(sum$chi,1) #biallelic, df=1
hwe<-sum[sum$chi.result > 0.05,]

#prune to keep only those found in most pops
sum.prune<-hwe[hwe$Pop=="ADULT" | hwe$Pop=="JUVIE" | 
	hwe$Pop == "POP",]
sum.sum<-tapply(sum.prune$N,sum.prune$Locus,sum)
sum.sum<-sum.sum[as.numeric(sum.sum) > 393]
sum.prune<-hwe[hwe$Locus %in% names(sum.sum),]

##prune based on representation in the groups
sum.list<-split(sum.prune, sum.prune$Pop)
adt.n<-sum.list$ADULT[sum.list$ADULT$N > 216 & !is.na(sum.list$ADULT$Hs),]
juv.n<-sum.list$JUVIE[sum.list$JUVIE$N > 157 & !is.na(sum.list$JUVIE$Hs),]
fem.n<-sum.list$FEM[sum.list$FEM$N > 57& !is.na(sum.list$FEM$Hs),]
mal.n<-sum.list$MAL[sum.list$MAL$N>159& !is.na(sum.list$MAL$Hs),]
mom.n<-sum.list$MOM[sum.list$MOM$N>133& !is.na(sum.list$MOM$Hs),]
prg.n<-sum.list$PREGGER[sum.list$PREGGER$N>159& !is.na(sum.list$PREGGER$Hs),]
non.n<-sum.list$NONPREG[sum.list$NONPREG$N>14& !is.na(sum.list$NONPREG$Hs),]

#Now prune based on allele frequencies
adt.n<-adt.n[adt.n$Allele1Freq > 0.05 & adt.n$Allele1Freq < 0.95,]
juv.n<-juv.n[juv.n$Allele1Freq > 0.05 & juv.n$Allele1Freq < 0.95,]
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
mal.n<-mal.n[mal.n$Allele1Freq > 0.05 & mal.n$Allele1Freq < 0.95,]
mom.n<-mom.n[mom.n$Allele1Freq > 0.05 & mom.n$Allele1Freq < 0.95,]
prg.n<-prg.n[prg.n$Allele1Freq > 0.05 & prg.n$Allele1Freq < 0.95,]
non.n<-non.n[non.n$Allele1Freq > 0.05 & non.n$Allele1Freq < 0.95,]

sum.keep<-sum[sum$Locus %in% adt.n$Locus & sum$Locus %in% juv.n$Locus
	& sum$Locus %in% fem.n$Locus & sum$Locus %in% mal.n$Locus
	& sum$Locus %in% mom.n$Locus,]#31027



#####THIS IS WHERE I BOOTSTRAP IT#####
#split it by chromosomes
sum.keep$Chrom<-factor(sum.keep$Chrom)
sum.split<-split(sum.keep,sample())

bootstrap.fsts<-function(x){
	nloci<-length(levels(as.factor(x$Locus)))
	print(paste(nloci," loci on chrom ",levels(factor(x$Chrom))))
	loci.idx<-sample(1:nloci,nloci,replace=T)
	loci<-levels(as.factor(x$Locus))
	loci.keep<-loci[loci.idx]
	d<-data.frame()
	for(i in 1:nloci){
		k<-x[x$Locus == loci.keep[i],]
		adt.juv.ht<-2*
			(((k[k$Pop=="ADULT","N"]*k[k$Pop=="ADULT","Allele1Freq"])+
			(k[k$Pop=="JUVIE","N"]*k[k$Pop=="JUVIE","Allele1Freq"]))/
			(k[k$Pop=="ADULT","N"]+k[k$Pop=="JUVIE","N"]))*
			(((k[k$Pop=="ADULT","N"]*k[k$Pop=="ADULT","Allele2Freq"])+
			(k[k$Pop=="JUVIE","N"]*k[k$Pop=="JUVIE","Allele2Freq"]))/
			(k[k$Pop=="ADULT","N"]+k[k$Pop=="JUVIE","N"]))
		fem.mal.ht<-2*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele1Freq"])+
			(k[k$Pop=="MAL","N"]*k[k$Pop=="MAL","Allele1Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MAL","N"]))*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele2Freq"])+
			(k[k$Pop=="MAL","N"]*k[k$Pop=="MAL","Allele2Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MAL","N"]))
		fem.mom.ht<-2*(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele1Freq"])+
			(k[k$Pop=="MOM","N"]*k[k$Pop=="MOM","Allele1Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MOM","N"]))*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele2Freq"])+
			(k[k$Pop=="MOM","N"]*k[k$Pop=="MOM","Allele2Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MOM","N"]))
		adt.juv.fst<-(adt.juv.ht-
			((k[k$Pop=="ADULT","Hs"]+k[k$Pop=="JUVIE","Hs"])/2))/adt.juv.ht
		fem.mal.fst<-(fem.mal.ht-
			((k[k$Pop=="FEM","Hs"]+k[k$Pop=="MAL","Hs"])/2))/fem.mal.ht
		fem.mom.fst<-(fem.mom.ht-
			((k[k$Pop=="FEM","Hs"]+k[k$Pop=="MOM","Hs"])/2))/fem.mom.ht
		d<-rbind(d,cbind(as.character(levels(factor(k$Chrom))),
			as.numeric(levels(factor(k$Pos))),
			as.character(levels(factor(k$Locus))),adt.juv.fst,
			fem.mal.fst,fem.mom.fst))
	}
	colnames(d)<-c("Chrom","Pos","Locus","AdultJuvieFst","FemMalFst",
		"FemMomFst")
	return(d)
}

boot.fsts<-lapply(sum.split,bootstrap.fsts)

#####bootstrap all fsts

bootstrap.all.fsts<-function(x){
	nloci<-length(levels(as.factor(x$Locus)))
	loci.idx<-sample(1:nloci,nloci,replace=T)
	loci<-levels(as.factor(x$Locus))
	loci.keep<-loci[loci.idx]
	d<-data.frame()
	for(i in 1:nloci){
		k<-x[x$Locus == loci.keep[i],]
		adt.juv.ht<-2*
			(((k[k$Pop=="ADULT","N"]*k[k$Pop=="ADULT","Allele1Freq"])+
			(k[k$Pop=="JUVIE","N"]*k[k$Pop=="JUVIE","Allele1Freq"]))/
			(k[k$Pop=="ADULT","N"]+k[k$Pop=="JUVIE","N"]))*
			(((k[k$Pop=="ADULT","N"]*k[k$Pop=="ADULT","Allele2Freq"])+
			(k[k$Pop=="JUVIE","N"]*k[k$Pop=="JUVIE","Allele2Freq"]))/
			(k[k$Pop=="ADULT","N"]+k[k$Pop=="JUVIE","N"]))
		fem.mal.ht<-2*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele1Freq"])+
			(k[k$Pop=="MAL","N"]*k[k$Pop=="MAL","Allele1Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MAL","N"]))*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele2Freq"])+
			(k[k$Pop=="MAL","N"]*k[k$Pop=="MAL","Allele2Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MAL","N"]))
		fem.mom.ht<-2*(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele1Freq"])+
			(k[k$Pop=="MOM","N"]*k[k$Pop=="MOM","Allele1Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MOM","N"]))*
			(((k[k$Pop=="FEM","N"]*k[k$Pop=="FEM","Allele2Freq"])+
			(k[k$Pop=="MOM","N"]*k[k$Pop=="MOM","Allele2Freq"]))/
			(k[k$Pop=="FEM","N"]+k[k$Pop=="MOM","N"]))
		adt.juv.fst<-(adt.juv.ht-
			((k[k$Pop=="ADULT","Hs"]+k[k$Pop=="JUVIE","Hs"])/2))/adt.juv.ht
		fem.mal.fst<-(fem.mal.ht-
			((k[k$Pop=="FEM","Hs"]+k[k$Pop=="MAL","Hs"])/2))/fem.mal.ht
		fem.mom.fst<-(fem.mom.ht-
			((k[k$Pop=="FEM","Hs"]+k[k$Pop=="MOM","Hs"])/2))/fem.mom.ht
		d<-rbind(d,cbind(as.character(levels(factor(k$Chrom))),
			as.numeric(levels(factor(k$Pos))),
			as.character(levels(factor(k$Locus))),adt.juv.fst,
			fem.mal.fst,fem.mom.fst))
	}
	colnames(d)<-c("Chrom","Pos","Locus","AdultJuvieFst","FemMalFst",
		"FemMomFst")
	return(d)
}

all.fsts<-bootstrap.all.fsts(sum.keep)



