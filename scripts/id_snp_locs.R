setwd("~/Projects/scovelli_genome")

gff<-read.delim("annotated_scovelli_scaffolds_March2016.gff",stringsAsFactors = F,header=F)
gagp<-read.delim("SSC_genome.agp",stringsAsFactors = F,skip=3,header=F)
sagp<-read.delim("SSC_scaffolds.agp",stringsAsFactors = F,skip=3,header=F)
# ggff<-read.delim("ssc_122016_chromlevel.gff",stringsAsFactors = F,header=F)
# sgff<-read.delim("ssc_122016_scafflevel.gff",stringsAsFactors = F,header=F)
# gff<-rbind(ggff,sgff)
# gagp<-read.delim("ssc_122016_chromlevel.agp",stringsAsFactors = F,skip=3,header=F)
# sagp<-read.delim("ssc_122016_scafflevel.agp",stringsAsFactors = F,skip=3,header=F)
agp<-rbind(gagp,sagp)
drad<-read.table("~/Projects/SCA/results/biallelic/keep.vcf",header=F)
poly.loc<-drad[,1:2]#use keep.vcf

#match polymorphic loci to genome
poly.agp<-function(poly, agp,gff){
  #poly is a row of a df with chrom, pos as rows 1,2
  #agp is an agp file.
  chrom<-as.factor(unlist(poly[1]))
  region<-agp[agp[,1] %in% chrom & as.numeric(agp[,2])<=as.numeric(poly[2]) & as.numeric(agp[,3]) >=as.numeric(poly[2]),]
  reg<-data.frame(Chrom=region[,1],Pos=unlist(poly[2]),Chrom.start=as.numeric(region[,2]),Chrom.end=as.numeric(region[,3]),
                        Component=region[,6],comp.start=as.numeric(region[,7]),comp.end=as.numeric(region[,8]),stringsAsFactors=F)
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
type[type=="CDS,contig,exon,gene,mRNA"]<-"exon"
type[type=="CDS,contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type[type=="contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type[type=="contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type[type=="contig,five_prime_UTR"]<-"five_prime_UTR"
type[type=="contig,exon,five_prime_UTR,gene,mRNA,three_prime_UTR"]<-"gene"
type[type=="contig,gene,mRNA"]<-"gene"
type[is.na(type)]<-"gap"
barplot(summary(factor(type)))

fm.sig<-read.delim("~/Projects/SCA/results/biallelic/fm.sig.txt")
fm.sig$comploc<-gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\1\\.\\3",fm.sig$SNP)
fm.loc<-data.frame(Chrom=gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\1\\",fm.sig$SNP),
                   Pos=gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\3",fm.sig$SNP),stringsAsFactors = FALSE)
fm.annotate<-do.call("rbind",apply(fm.loc,1,poly.agp,agp=agp,gff=gff))
type.fm<-fm.annotate[,c("Chrom","Pos","type")]
type.fm$type[type.fm$type=="CDS,contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.fm$type[type.fm$type=="CDS,contig,exon,gene,mRNA"]<-"gene"
type.fm$type[type.fm$type=="CDS,contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.fm$type[type.fm$type=="contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.fm$type[type.fm$type=="contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.fm$type[type.fm$type=="contig,five_prime_UTR"]<-"five_prime_UTR"
type.fm$type[type.fm$type=="contig,exon,five_prime_UTR,gene,mRNA,three_prime_UTR"]<-"gene"
type.fm$type[type.fm$type=="contig,gene,mRNA"]<-"gene"
type.fm$type[is.na(type.fm$type)]<-"gap"
barplot(summary(factor(type.fm$type)))

mo.sig<-read.delim("~/Projects/SCA/results/biallelic/mo.sig.txt")
mo.sig$comploc<-gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\1\\.\\3",mo.sig$SNP)
mo.loc<-data.frame(Chrom=gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\1\\",mo.sig$SNP),
                   Pos=gsub("(\\w+\\d+)\\.(\\d+)\\.(\\d+)","\\3",mo.sig$SNP),stringsAsFactors = FALSE)
mo.annotate<-do.call("rbind",apply(mo.loc,1,poly.agp,agp=agp,gff=gff))
type.mo<-mo.annotate[,c("Chrom","Pos","type")]
type.mo$type[type.mo$type=="CDS,contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.mo$type[type.mo$type=="CDS,contig,exon,gene,mRNA"]<-"gene"
type.mo$type[type.mo$type=="CDS,contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.mo$type[type.mo$type=="contig,exon,five_prime_UTR,gene,mRNA"]<-"five_prime_UTR"
type.mo$type[type.mo$type=="contig,exon,gene,mRNA,three_prime_UTR"]<-"three_prime_UTR"
type.mo$type[type.mo$type=="contig,five_prime_UTR"]<-"five_prime_UTR"
type.mo$type[type.mo$type=="contig,exon,five_prime_UTR,gene,mRNA,three_prime_UTR"]<-"gene"
type.mo$type[type.mo$type=="contig,gene,mRNA"]<-"gene"
type.mo$type[is.na(type.mo$type)]<-"gap"
barplot(summary(factor(type.mo$type)))

bayes.test<-data.frame(Dataset=c(rep("SCA",(nrow(type.mo)+nrow(type.fm)))),
                       Species=c(rep("S.scovelli",(nrow(type.mo)+nrow(type.fm)))),
                       SelectionType=c(rep("SexualSelection",(nrow(type.mo))),
                                       rep("MFViabiliatySelection",nrow(type.fm))),
                       Habitat=c(rep("Marine",(nrow(type.mo)+nrow(type.fm)))),
                       Chrom=c(type.mo$Chrom,type.fm$Chrom),Pos=c(type.mo$Pos,type.fm$Pos),
                       type=c(type.mo$type,type.fm$type),stringsAsFactors=FALSE)