#Author: Sarah P. Flanagan
#Date: 9 April 2016
#Purpose: to analyze the outlier blast2go results.

rm(list=ls())
library(ggplot2)
setwd("E:/ubuntushare/SCA/results/biallelic_outliers/rad_region/blast2go")

go.plot<-function(file.list, file.name,analysis.list=NULL,pdf=FALSE){
	dat<-read.table(file.list[1],skip=1,sep='\t')
	dat<-dat[,1:2]
	if(!is.null(analysis.list)){
		analysis.names<-analysis.list
	} else {
		analysis.names<-gsub("(\\w+)_\\w+.*","\\1",file.list)
	}
	colnames(dat)<-c("GO","Number")
	dat$Analysis<-analysis.names[1]
	for(i in 2:length(file.list)){
		d<-read.table(file.list[i],skip=1,sep='\t')
		d<-d[,1:2]
		colnames(d)<-c("GO","Number")
		d$Analysis<-analysis.names[i]
		dat<-rbind(dat,d)
	}

	for(i in 1:length(levels(dat$GO))){
		t<-dat[dat$GO %in% levels(dat$GO)[i],]
		if(nrow(t)<length(analysis.names)){
			a.new<-analysis.names[!(analysis.names %in% t$Analysis)]
			g.new<-rep(t$GO[1],length(a.new))
			n.new<-rep(0,length(a.new))
			t.add<-data.frame(GO=g.new,Number=n.new,Analysis=a.new)
			dat<-rbind(dat,t.add)
		}
	}
	dat<-dat[order(dat$GO),]
	library(ggplot2)
	if(pdf==T){ 
		file.name<-paste(file.name,'.pdf',sep="") 
		pdf(file.name,height=9,width=9)
	}else {	
		file.name<-paste(file.name,'.jpeg',sep="") 
		jpeg(file.name,height=9,width=9,units="in",res=300)
	}
	par(mar=c(2,2,2,2),oma=c(2,2,2,2),cex=2,lwd=1.3)
	p<-ggplot(dat,aes(factor(GO),Number,fill = factor(Analysis))) + 
		geom_bar(stat="identity",position="dodge") + 
		scale_fill_brewer(palette="Set1",name="Analysis") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		coord_flip() +
		xlab("Gene Ontology") + ylab("Number")
	print(p)
	dev.off()
	return(dat)
}	
unique<-c("aj_","fm_","mo_","lrt_mo","lrt_fm")
bio.files<-list.files(pattern="biol.txt")
bio2.files<-list.files(pattern="biol2.txt")
cell.files<-list.files(pattern="cell.txt")
cell2.files<-list.files(pattern="cell2.txt")
mol.files<-list.files(pattern="mol.txt")
mol2.files<-list.files(pattern="mol2.txt")

bio.unique<-bio.files[sub("(\\w{2}_)\\w+.*","\\1",bio.files) %in% unique]
bio2.unique<-bio2.files[sub("(\\w{2}_)\\w+.*","\\1",bio2.files) %in% unique]
cell.unique<-cell.files[sub("(\\w{2}_)\\w+.*","\\1",cell.files) %in% unique]
cell2.unique<-cell2.files[sub("(\\w{2}_)\\w+.*","\\1",cell2.files) %in% unique]
mol.unique<-mol.files[sub("(\\w{2}_)\\w+.*","\\1",mol.files) %in% unique]
mol2.unique<-mol2.files[sub("(\\w{2}_)\\w+.*","\\1",mol2.files) %in% unique]

comparisons<-c("ajfm_","fmmo_","shared_")
bio.comp<-bio.files[sub("(\\w{2}_)\\w+.*","\\1",bio.files) %in% comparisons]
bio2.comp<-bio2.files[sub("(\\w{2}_)\\w+.*","\\1",bio2.files) %in% comparisons]
cell.comp<-cell.files[sub("(\\w{2}_)\\w+.*","\\1",cell.files) %in% comparisons]
cell2.comp<-cell2.files[sub("(\\w{2}_)\\w+.*","\\1",cell2.files) %in% comparisons]
mol.comp<-mol.files[sub("(\\w{2}_)\\w+.*","\\1",mol.files) %in% comparisons]
mol2.comp<-mol2.files[sub("(\\w{2}_)\\w+.*","\\1",mol2.files) %in% comparisons]


analysis.names<-c("Adult-Offspring","Males-Females","Mothers-Females","Shared")
bio.dat<-go.plot(bio.comp,"../../Biology",analysis.names,pdf=TRUE)
bio2.dat<-go.plot(bio2.comp,"Biology2",analysis.names)
cell.dat<-go.plot(cell.comp,"Cell",analysis.names)
cell2.dat<-go.plot(cell2.comp,"Cell2",analysis.names)
mol.dat<-go.plot(mol.comp,"Molecular",analysis.names)
mol2.dat<-go.plot(mol2.comp,"Molecular2",analysis.names)

#Analyze which ones are alone.
bio.split<-split(bio.dat,bio.dat$GO)
bio.analysis<-data.frame(Status=unlist(lapply(bio.split,function(x){
		out<-"NotShared"
		if(nrow(x[x$Number>0,])==1){
			out<-x$Analysis[x$Number>0]
		}
		return(out) })))
