#Try doing gwsca_biallelic_vcf here.

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")
calc.fst<-function(x,y){ #x is one row of a vcf file
	if(x[1]==y[1] & x[2]==y[2]){
	chr=x[1]
	pos=x[2]	
	x<-unlist(x[3:length(x)])
	x<-factor(x[x!="./."])
	y<-unlist(y[3:length(y)])
	y<-factor(y[y!="./."])
	xy<-c(as.character(x),as.character(y))
	gf.x<-table(x)/sum(table(x))
	af.x<-table(unlist(lapply(as.character(x),strsplit,"/")))/
		sum(table(unlist(lapply(as.character(x),strsplit,"/"))))
	hs.x<-2*af.x[1]*af.x[2]
	gf.y<-table(y)/sum(table(y))
	af.y<-table(unlist(lapply(as.character(y),strsplit,"/")))/
		sum(table(unlist(lapply(as.character(y),strsplit,"/"))))
	hs.y<-2*af.y[1]*af.y[2]
	if(length(af.x)>1 & length(af.y) > 1){
	af<-table(unlist(lapply(as.character(xy),strsplit,"/")))/
		sum(table(unlist(lapply(as.character(xy),strsplit,"/"))))
	hs<-(hs.x*length(x)+hs.y*length(y))/(length(x)+length(y))
	ht<-2*af[1]*af[2]
	fst<-(ht-hs)/ht
	return(data.frame(Chr=chr,Pos=pos,MajAF1=max(af.x),
		MajAF2=max(af.y), N1=length(x),N2=length(y),Hs=hs,Ht=ht,Fst=fst))
	} else {
		return(data.frame(Chr=chr,Pos=pos,MajAF1=NA,
		MajAF2=NA, N1=length(x),N2=length(y),Hs=NA,Ht=NA,Fst=NA))
	}
} else {
	print("Sort your vcfs! The two do not match up.")
}}
prune.vcf<-function(vcf.df,cov.per){
	keep<-apply(vcf.df,1,function(x) length(x[x!="./."]))
	names(keep)<-seq(1,length(keep))
	keepnum<-cov.per*(ncol(vcf.df)-2)
	vcf.df<-vcf.df[names(keep[keep > keepnum]),]
	return(vcf.df)
}

vcf.orig<-read.delim("biallelic.gt.vcf")
info<-read.delim("ind_info_vcf.txt",col.names=c("name","ID","sex","age","status"))

#prune for overall coverage
vcf<-prune.vcf(vcf.orig,0.75)

adt<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$age=="ADULT",]$name])
off<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$age=="JUVIE",]$name])
fem<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$sex=="FEM",]$name])
mal<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$sex=="MAL",]$name])
mom<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$status=="MOM",]$name])

#calculate fsts, then prune
adt<-adt[order(c(adt$CHROM,adt$POS)),]
off<-off[order(c(off$CHROM,off$POS)),]
mal<-mal[order(c(mal$CHROM,mal$POS)),]
fem<-fem[order(c(fem$CHROM,fem$POS)),]
mom<-mom[order(c(mom$CHROM,mom$POS)),]

ao.fsts<-data.frame()
for(i in 1:nrow(vcf)){
	ao.fsts<-rbind(ao.fsts,calc.fst(adt[i,],off[i,]))
}

fm.fsts<-data.frame()
for(i in 1:nrow(vcf)){
	fm.fsts<-rbind(fm.fsts,calc.fst(mal[i,],fem[i,]))
}

md.fsts<-data.frame()
for(i in 1:nrow(vcf)){
	md.fsts<-rbind(md.fsts,calc.fst(mom[i,],fem[i,]))
}

#Remove any NAs
ao.prune<-ao.fsts[!is.na(ao.fsts$Fst),]
fm.prune<-fm.fsts[!is.na(fm.fsts$Fst),]
md.prune<-md.fsts[!is.na(md.fsts$Fst),]
#Prune for per-group Coverage
ao.prune<-ao.prune[ao.prune$N1 >= 0.5*(ncol(adt)-2) &
	ao.prune$N2 >= 0.5*(ncol(off)-2),]
fm.prune<-fm.prune[fm.prune$N1 >= 0.5*(ncol(mal)-2) &
	fm.prune$N2 >= 0.5*(ncol(fem)-2),]
md.prune<-md.prune[md.prune$N1 >= 0.5*(ncol(mom)-2) &
	md.prune$N2 >= 0.5*(ncol(fem)-2),]

#Prune for Allele Frequency
ao.prune<-ao.prune[ao.prune$MajAF1 > 0.05 & ao.prune$MajAF1 < 0.95 &
	ao.prune$MajAF2 > 0.05 & ao.prune$MajAF2 < 0.95,]
fm.prune<-fm.prune[fm.prune$MajAF1 > 0.05 & fm.prune$MajAF1 < 0.95 &
	fm.prune$MajAF2 > 0.05 & fm.prune$MajAF2 < 0.95,]
md.prune<-md.prune[md.prune$MajAF1 > 0.05 & md.prune$MajAF1 < 0.95 &
	md.prune$MajAF2 > 0.05 & md.prune$MajAF2 < 0.95,]


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

#plot with the model CIs
png("fst.biallelic.pruned.Rcalc.model.png",height=300,width=300,units="mm",res=300)
par(mfrow=c(3,1),oma=c(1,1,0,0),mar=c(0,1,1,0),mgp=c(3,0.5,0), cex=1.5)
plot.fsts(ao.prune, ci.dat=aj.null,fst.name="Fst", chrom.name="CHROM"
	, axis.size=0.75, bp.name="POS")
legend("top","Adult-Juvenile", cex=0.75,bty="n")
plot.fsts(fm.prune, ci.dat=mf.null,fst.name="Fst", chrom.name="CHROM"
	, axis.size=0.75,bp.name="POS")
legend("top","Male-Female", cex=0.75,bty="n")
plot.fsts(md.prune, ci.dat=mo.null,fst.name="Fst", chrom.name="CHROM"
	, axis.size=0.75,bp.name="POS")
legend("top","Mothers-Females", cex=0.75,bty="n")
mtext("Genomic Location", 1, outer=T, cex=1)
mtext("Fst", 2, outer=T, cex=1)
dev.off()

