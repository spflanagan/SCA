#Author: Sarah P. Flanagan
#Last Updated: 11 May 2016
#Purpose: Calculate the error rate based on consistency of genotyping the same
#individuals.

setwd("E:/ubuntushare/SCA/results/biallelic")
vcf<-read.table("biallelic_merge.vcf",header=T)
duped<-c("sample_OFF016.1_align","sample_OFF016_align",
	"sample_OFF027.1_align","sample_OFF027_align",
	"sample_OFF032.1_align","sample_OFF032_align",
	"sample_PRM177.1_align","sample_PRM177_align")
vcf.compare<-cbind(vcf[,1:9],vcf[,duped])

compare.two.inds<-function(vcf,column1,column2){
	all.gen1<-nrow(vcf[vcf[,column1]!= "./.",])
	all.gen2<-nrow(vcf[vcf[,column2]!= "./.",])
	same.gen<-nrow(vcf[vcf[,column1] == vcf[,column2] & 
		vcf[,column1] != "./.",])
	same.missing<-nrow(vcf[vcf[,column1] == vcf[,column2] & 
		vcf[,column1] == "./.",])
	missing.one<-nrow(vcf[vcf[,column1] != vcf[,column2] & 
		vcf[,column1] == "./.",])
	missing.two<-nrow(vcf[vcf[,column1] != vcf[,column2] & 
		vcf[,column2] == "./.",])
	diff.gen<-nrow(vcf[vcf[,column1] != vcf[,column2] & 
		vcf[,column1] != "./." & vcf[,column2] != "./.",])
	summary<-data.frame(Same.Genotype=same.gen,Same.Missing=same.missing,
		Different.Genotypes=diff.gen,
		Missing.In.One=missing.one,Missing.In.Other=missing.two,
		Total.Missing.In.One=missing.one+missing.two,
		Num.Genotyped.One=all.gen1, Num.Genotyped.Other=all.gen2,
		Number.Loci=nrow(vcf))
	rownames(summary)<-colnames(vcf)[column2]
	return(summary)
}

error.dat<-data.frame()
for(i in seq(10,ncol(vcf.compare),2)){
	sumi<-compare.two.inds(vcf.compare,i,i+1)
	error.dat<-rbind(error.dat,sumi)
}
error.dat<-error.dat/error.dat$Number.Loci
error.dat$Number.Loci<-nrow(vcf.compare)
write.table(error.dat, "ErrorRates_allLoci.txt",col.names=T,row.names=T,quote=F)

#What about those loci that are in the analyses?
prg.n.loc<-read.table("LociInPregMales.txt",header=T)
juv.n.loc<-read.table("LociInOffspring.txt",header=T)
vcf.compare$Locus<-paste(vcf.compare$CHROM,vcf.compare$ID,vcf.compare$POS,sep=".")
prg.n.loc$Locus<-paste(prg.n.loc$Chrom,prg.n.loc$LocID,prg.n.loc$Pos,sep=".")
juv.n.loc$Locus<-paste(juv.n.loc$Chrom,juv.n.loc$LocID,juv.n.loc$Pos,sep=".")
overlapping<-prg.n.loc[prg.n.loc$Locus %in% juv.n.loc$Locus,]
vcf.juv<-vcf.compare[vcf.compare$Locus %in% juv.n.loc$Locus,]
vcf.prg<-vcf.compare[vcf.compare$Locus %in% prg.n.loc$Locus,]


prune.error.dat<-data.frame()
for(i in seq(10,(ncol(vcf.juv)-3),2)){
	sumi<-compare.two.inds(vcf.juv,i, i+1)
	prune.error.dat<-rbind(prune.error.dat,sumi)
}
sumi<-compare.two.inds(vcf.prg,16,17)
prune.error.dat<-rbind(prune.error.dat,sumi)
nloci<-prune.error.dat$Number.Loci

prune.error.dat<-prune.error.dat/prune.error.dat$Number.Loci
prune.error.dat$Number.Loci<-nloci
write.table(prune.error.dat, "ErrorRates_prunedLoci.txt",col.names=T,
	row.names=T,quote=F)

