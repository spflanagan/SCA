setwd("~/Projects/scovelli_genome")

gff<-read.delim("annotated_scovelli_scaffolds_March2016.gff",stringsAsFactors = F,header=F)
gagp<-read.delim("SSC_genome.agp",stringsAsFactors = F,skip=3,header=F)
sagp<-read.delim("SSC_scaffolds.agp",stringsAsFactors = F,skip=3,header=F)
agp<-rbind(gagp,sagp)
poly.loc<-drad[,c("#CHROM","POS")]#use keep.vcf

#match polymorphic loci to genome
poly.agp<-function(poly, agp,gff){
  #poly is a row of a df with chrom, pos as rows 1,2
  #agp is an agp file.
  chrom<-as.factor(unlist(poly[1]))
  region<-agp[agp[,1] %in% chrom & as.numeric(agp[,2])<=as.numeric(poly[2]) & as.numeric(agp[,3]) >=as.numeric(poly[2]),]
  reg<-data.frame(stringsAsFactors=F,
                  cbind(Chrom=region[,1],Pos=poly[2],Chrom.start=as.numeric(region[,2]),Chrom.end=as.numeric(region[,3]),
                        Component=region[,6],comp.start=as.numeric(region[,7]),comp.end=as.numeric(region[,8])))
  adj.pos<-(as.numeric(reg$Pos)-as.numeric(reg$Chrom.start))+as.numeric(reg$comp.start)-1
  if(nrow(gff[gff$V1 %in% reg$Component,])>0){
    chrom<-reg$Component
  }
  if(nrow(gff[gff$V1 %in% reg$Chrom,])>0){
    chrom<-reg$Chrom
  }
  g<-gff[gff$V1 %in% chrom & gff$V4 <= adj.pos & gff$V5 >= adj.pos,]
  type<-levels(as.factor(g$V3))
  t<-type[1]
  if(length(type) > 1){
    for(i in 2:length(type)){
      t<-paste(t,type[i],sep=",")
    }
  }
  reg$adj.pos<-adj.pos
  reg$type<-t
  return(reg)
}

drad.annotate<-do.call("rbind",apply(poly.loc,1,poly.agp,agp=agp,gff=gff))
type<-drad.annotate$type
type[type=="CDS,contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type[type=="CDS,contig,exon,gene,mRNA"]<-"gene"
type[type=="CDS,contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type[type=="contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type[type=="contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type[type=="contig,five_prime_UTR"]<-"five_prime_UTR"
type[type=="contig,exon,five_prime_UTR,gene,mRNA,three_prime_UTR"]<-"gene"
type[type=="contig,gene,mRNA"]<-"gene"
type[is.na(type)]<-"gap"
barplot(summary(factor(type)))

type.fm<-drad.annotate$type[paste(drad.annotate$Chrom,drad.annotate$Pos,sep=".") %in% fm.sig$comploc]
type.fm[type.fm=="CDS,contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.fm[type.fm=="CDS,contig,exon,gene,mRNA"]<-"gene"
type.fm[type.fm=="CDS,contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.fm[type.fm=="contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.fm[type.fm=="contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.fm[type.fm=="contig,five_prime_UTR"]<-"five_prime_UTR"
type.fm[type.fm=="contig,exon,five_prime_UTR,gene,mRNA,three_prime_UTR"]<-"gene"
type.fm[type.fm=="contig,gene,mRNA"]<-"gene"
type.fm[is.na(type.fm)]<-"gap"
barplot(summary(factor(type.fm)))
