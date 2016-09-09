#Author: Sarah P. Flanagan
#Date: 7 September 2016
#Purpose: to analyze the outlier blast2go results.

getSrcDirectory(function(x) {x})
rm(list=ls())
library(ggplot2)
#setwd("B:/ubuntushare/SCA/")
setwd("~/Projects/SCA/results/biallelic_outliers/rad_region/blast2go")

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
	dat$Number<-dat$Number/sum(dat$Number)
	for(i in 2:length(file.list)){
		d<-read.table(file.list[i],skip=1,sep='\t')
		d<-d[,1:2]
		colnames(d)<-c("GO","Number")
		d$Analysis<-analysis.names[i]
		d$Number<-d$Number/sum(d$Number)
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
		pdf(file.name,height=10,width=9)
	}else {	
		file.name<-paste(file.name,'.jpeg',sep="") 
		jpeg(file.name,height=10,width=9,units="in",res=300)
	}
	par(mar=c(2,2,2,2),oma=c(2,2,2,2),cex=2,lwd=1.3)
	p<-ggplot(dat,aes(factor(GO),Number,fill = factor(Analysis))) + 
		geom_bar(stat="identity",position="dodge") + 
		scale_fill_brewer(palette="Set1",name="Analysis") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		coord_flip() +
		xlab("Gene Ontology") + ylab("Proportion")
	print(p)
	dev.off()
	return(dat)
}	


bio.files<-list.files(pattern="biol.txt")
bio2.files<-list.files(pattern="biol2.txt")
bio3.files<-list.files(pattern="biol3.txt")
cell.files<-list.files(pattern="cell.txt")
cell2.files<-list.files(pattern="cell2.txt")
mol.files<-list.files(pattern="mol.txt")
mol2.files<-list.files(pattern="mol2.txt")

unique<-c("aj_","fm_","mo_","lrt_mo","lrt_fm","sharedall")
comparisons<-c("aj","fm","mo","lrt_mo","lrt_fm","sharedall")

bio.unique<-bio.files[sub("(\\w{2}_)\\w+.*","\\1",bio.files) %in% unique]
bio2.unique<-bio2.files[sub("(\\w{2}_)\\w+.*","\\1",bio2.files) %in% unique]
cell.unique<-cell.files[sub("(\\w{2}_)\\w+.*","\\1",cell.files) %in% unique]
cell2.unique<-cell2.files[sub("(\\w{2}_)\\w+.*","\\1",cell2.files) %in% unique]
mol.unique<-mol.files[sub("(\\w{2}_)\\w+.*","\\1",mol.files) %in% unique]
mol2.unique<-mol2.files[sub("(\\w{2}_)\\w+.*","\\1",mol2.files) %in% unique]
bio.comp<-bio.files[sub("(\\w+)_biol.txt","\\1",bio.files) %in% comparisons]
bio2.comp<-bio2.files[sub("(\\w+)_biol2.txt","\\1",bio2.files) %in% comparisons]
bio3.comp<-bio3.files[sub("(\\w+)_biol3.txt","\\1",bio3.files) %in% comparisons]
cell.comp<-cell.files[sub("(\\w+)_cell.txt","\\1",cell.files) %in% comparisons]
cell2.comp<-cell2.files[sub("(\\w+)_cell2.txt","\\1",cell2.files)%in% comparisons]
mol.comp<-mol.files[sub("(\\w+)_mol.txt","\\1",mol.files) %in% comparisons]
mol2.comp<-mol2.files[sub("(\\w+)_mol2.txt","\\1",mol2.files) %in% comparisons]

#get the total number of hits to add to the analysis names
n.bio2<-NULL
for(i in 1:length(bio2.comp))
{
  dat<-read.table(bio2.comp[i],skip=1,sep='\t')
  n.bio2[i]<-sum(dat$V2)
  names(n.bio2)[i]<-bio2.comp[i]
}
#aj_biol3.txt        fm_biol3.txt    lrt_fm_biol3.txt    lrt_mo_biol3.txt        mo_biol3.txt sharedall_biol3.txt 
#891                 665                 104                  19                 464                 205 
#aj_biol2.txt        fm_biol2.txt    lrt_fm_biol2.txt    lrt_mo_biol2.txt        mo_biol2.txt sharedall_biol2.txt 
#746                 527                  76                  13                 376                 175 
> 
analysis.names<-c("Fst Adult-Offspring (746)","Fst Males-Females (527)",
	"LRT Males-Females (76)","LRT Mothers-Adults (13)","Fst Mothers-Females (376)",
	"Fst Shared (175)")

bio.dat<-go.plot(bio.comp,"Biology",analysis.names)
bio2.dat<-go.plot(bio2.comp,"Biology2",analysis.names)
bio3.dat<-go.plot(bio3.comp,"Biology3",analysis.names)
bio3.dat<-go.plot(bio3.comp,"Biology3",analysis.names,pdf=T)
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

#Get the blast2go tables
table.files<-list.files(pattern="blast2go")
table.files<-table.files[sub("(\\w+)_blast2go.txt","\\1",table.files) %in% comparisons]
tables<-lapply(table.files,read.delim)
names(tables)<-table.files
info<-data.frame(Comparison=character(),PropShared=numeric(),PropUnique=numeric(),
                 PropNA=numeric(), NumGenes=numeric(),stringsAsFactors=F)
for(i in 1:length(tables)){
  gene1<-tables[[(i)]]$Description
  na<-length(gene1[gene1=="---NA---"])
  gene1<-gene1[gene1!="---NA---"]
  num.genes<-as.numeric(length(gene1[!duplicated(gene1)]))
  shared.genes<-NULL
  for(j in i:length(tables)){
    if(i!=j)
    {
      gene2<-tables[[j]]$Description
      gene2<-gene2[gene2!="---NA---"]
      shared.genes<-gene1[gene1 %in% gene2]
      shared.genes<-shared.genes[!duplicated(shared.genes)]
    }
  }
  info[i,]<-c(names(tables)[i],as.numeric(length(shared.genes))/num.genes,
              (num.genes-as.numeric(length(shared.genes)))/num.genes,na/num.genes,num.genes)
}

dat<-data.frame(Description=character(),Number=numeric(),Analysis=character(),stringsAsFactors = F)
x<-1
for(i in 1:nrow(info)){
  dat[x+0,]<-c(colnames(info)[2],info[i,2],info[i,1])
  dat[x+1,]<-c(colnames(info)[3],info[i,3],info[i,1])
  dat[x+2,]<-c(colnames(info)[4],info[i,4],info[i,1])
  x<-x+3
}
jpeg("SharedUniqueGenes",height=10,width=9,units="in",res=300)
par(mar=c(2,2,2,2),oma=c(2,2,2,2),cex=2,lwd=1.3)
p<-ggplot(dat,aes(factor(Description),Number,fill = factor(Analysis))) + 
  geom_bar(stat="identity",position="dodge") + 
  scale_fill_brewer(palette="Set1",name="Analysis") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() +
  xlab("") + ylab("Proportion")
print(p)
dev.off()


table<-read.delim(table.files[1])
find.bio.desc<-function(x){
  if(x!=""){
    t<-unlist(strsplit(as.character(x),";"))
    t<-unlist(lapply(t,gsub,pattern=" (\\w:)",replacement="\\1"))
    #find and keep only those with P
    k<-lapply(t,grep,pattern="P:")
    tk<-t[k==1]
    tk<-tk[!is.na(tk)]
    p<-unlist(lapply(tk, gsub,pattern="P:(\\w+)",replacement="\\1"))
  } else{
    p<-0
  }
  return(p)
}
all.ps<-unlist(lapply(table$GO.Names.list, find.bio.desc))
freq.all.ps<-table(all.ps)
ps<-lapply(table$GO.Names.list, find.bio.desc)
write.table("SeqName\tGOterms","BiolGOTerms.txt",col.names=F,quote=F)
for(i in 1:length(ps)){
  write.table(cbind(table$SeqName,"BiolGOTerms.txt",ps[i]),"BiolGOTerms.txt",append=T,col.names=F,quote=F,row.names=F,sep='\t')
  
}

