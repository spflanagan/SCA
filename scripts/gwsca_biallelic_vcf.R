#Try doing gwsca_biallelic_vcf here.

source("E:/ubuntushare/SCA/scripts/plotting_functions.R")
setwd("E:/ubuntushare/SCA/results/biallelic")
vcf<-read.delim("biallelic.gt.vcf")
info<-read.delim("ind_info_vcf.txt",col.names=c("name","ID","sex","age","status"))

adults<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$age=="ADULT",]$name])
off<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$age=="JUVIE",]$name])
fem<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$sex=="FEM",]$name])
mal<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$sex=="MAL",]$name])
mom<-cbind(vcf[,1:2],vcf[,colnames(vcf) %in% info[info$status=="MOM",]$name])

calc.fst<-function(x,y){ #x is one row of a vcf file
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
	af<-table(unlist(lapply(as.character(xy),strsplit,"/")))/
		sum(table(unlist(lapply(as.character(xy),strsplit,"/"))))
	hs<-(hs.x*length(x)+hs.y*length(y))/(length(x)+length(y))
	ht<-2*af[1]*af[2]
	fst<-(ht-hs)/ht
	return(data.frame(Hs=hs,Ht=ht,Fst=fst))
}

adults<-adults[order(c(adults$CHROM,adults$POS)),]
off<-off[order(c(off$CHROM,off$POS)),]

fsts<-data.frame()
for(i in 1:nrow(vcf)){
	fsts<-rbind(fsts,calc.fst(adults[i,],off[i,]))
}

mal<-mal[order(c(mal$CHROM,mal$POS)),]
fem<-fem[order(c(fem$CHROM,fem$POS)),]
fm.fsts<-data.frame()
for(i in 1:nrow(vcf)){
	fm.fsts<-rbind(fm.fsts,calc.fst(mal[i,],fem[i,]))
}

mom<-adults[order(c(mom$CHROM,mom$POS)),]
mom.fsts<-data.frame()
for(i in 1:nrow(vcf)){
	mom.fsts<-rbind(mom.fsts,calc.fst(fem[i,],mom[i,]))
}


