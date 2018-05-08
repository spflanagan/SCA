#Author: Sarah P. Flanagan
#Last Updated: 3 June 2016
#Date: 15 February 2016
#Purpose: Analyze the output from CERVUS when it was run with various numbers
#of loci

rm(list=ls())
#setwd("E:/ubuntushare/SCA/results/")
setwd("GitHub/SCA/results")

##############################################################################
#####CERVUS ANALYSIS
##############################################################################

##### FUNCTIONS #####
cervus_analysis<-function(stats.files,pattern,match_name="Mother given known father:"){
  stats<-data.frame(NumLoci=numeric(),ConfidenceLevel=numeric(), Delta=numeric(),
  	NumAssignments=numeric(), AssignmentRate=numeric(),
  	MarkerType=character(),stringsAsFactors=FALSE)
  for(i in 1: length(stats.files)){
  	dat<- readLines(stats.files[i])
  	num.loci<-as.numeric(gsub(pattern,"\\1",stats.files[i]))
  	dir<-gsub("(.*)/.*\\d+_\\d+.*.txt","\\1",stats.files[i])
  	start<-match(match_name, dat)
  	#pull out info for strict only
  	info<-unlist(strsplit(dat[(start+4)],"\\s+"))
  	info<-unlist(strsplit(gsub("[()%]","",info),"[[:space:]]"))
  	stats[i,]<-cbind(as.numeric(num.loci),as.numeric(info[2]),
  		as.numeric(info[3]),as.numeric(info[4]), as.numeric(info[6]),dir)
  	rownames(stats)[i]<-as.character(stats.files[i])
  }
  return(stats)
}

plot_delta<-function(stats,cols,borders,leg.loc){
  plot(c(0.5,16.5),c(min(as.numeric(stats$Delta)),max(as.numeric(stats$Delta))),type='n',axes=FALSE,xlab="",ylab="")
  abline(h=0,lty=2,col="darkgrey")
  boxplot(as.numeric(stats$Delta)~stats$MarkerType*as.numeric(stats$NumLoci),
          col=cols,notch=FALSE,add=TRUE,
          border=borders,las=1,
          ylab="",xlab="",xaxt='n')
  axis(1,at=seq(1.5,length(sort(as.numeric(unique(stats$NumLoci))))*2+0.5,2),
       labels=sort(as.numeric(unique(stats$NumLoci))),las=2)
  mtext(expression(Critical~Delta),2,outer=F, line=2)
  
  legend(leg.loc,c("Haplotypes","SNPs"),pt.bg=cols,
         bty='n',pch=22,col=borders)
}
plot_assignmentRate<-function(stats,cols,haps.name,snps.name,r,col2plot="AssignmentRate",
                              ylab="Assignment Rate (%)",loc.name="NumLoci",y.max=NULL){
  if(is.null(y.max)){
    y.max<-round((max(as.numeric(stats[stats$MarkerType==haps.name,col2plot]),
             as.numeric(stats[stats$MarkerType==snps.name,col2plot]))+10)/10)*10
  }
  plot(as.numeric(stats[stats$MarkerType==haps.name,loc.name])-10, 
       as.numeric(stats[stats$MarkerType==haps.name,col2plot]),
       xaxt='n',las=1,col=cols[1],ylab="",xlab="", pch=19,ylim=c(0,y.max),
       xlim=c(min(as.numeric(stats$NumLoci)),max(as.numeric(stats$NumLoci))))
  points(as.numeric(stats[stats$MarkerType==snps.name,loc.name])+10, 
         as.numeric(stats[stats$MarkerType==snps.name,col2plot]),
         col=cols[2],pch=15)
  axis(1, at=sort(as.numeric(unique(stats$NumLoci))),las=2)
  pts<-mapply(function(means,jig,cols){
    pts<-data.frame(x=as.numeric(names(means)),y=means)
    pts$xmin<-pts$x+jig-35
    pts$xmax<-pts$x+jig+35
    apply(pts,1,function(pt,color){
      lines(x=c(pt["xmin"],pt["xmax"]),y=rep(pt["y"],2),lwd=2,col=color)
    },color=cols)
    return(data.frame(pts))
  },means=r[c(haps.name,snps.name)],jig=c(-10,10),cols=cols)
  mtext(ylab,2,outer=F,line=2)
}

cervus_patcorr<-function(csv.files,file.patt,name.patt,off.keep=NULL){
  prop<-do.call(rbind,lapply(csv.files,function(filename,name.patt,file.patt){
    assigns<-read.csv(filename,row.names = NULL)
    colnames(assigns)<-colnames(assigns)[-1] #fix the column names
    assigns<-assigns[assigns$Offspring.ID %in% off.keep,]
    #delta.thresh<-pats[rownames(pats) == gsub("csv","txt",filename),"Delta"]
    sig.assigns<-assigns[assigns$Pair.confidence=="*",]
    if(length(grep("_align",sig.assigns$Offspring.ID))>0){#remove any tricky things
      sig.assigns$Offspring.ID<-gsub("(.*)_align","\\1",sig.assigns$Offspring.ID)
      sig.assigns$Candidate.father.ID<-gsub("(.*)_align","\\1",sig.assigns$Candidate.father.ID)
    }
    sig.assigns$offn<-gsub(paste(name.patt,"(-\\d.*)?",sep=""),"\\1",sig.assigns$Offspring.ID)
    sig.assigns$dadn<-gsub(paste(name.patt,"(.*)",sep=""),"\\1\\2",sig.assigns$Candidate.father.ID)
    sig.assigns$dadn[grep("086-23",sig.assigns$dadn)]<-"08623"
    sig.assigns$dadn[grep("086R",sig.assigns$dadn)]<-"086"
    prop.corr<-nrow(sig.assigns[sig.assigns$dadn==sig.assigns$offn,])/nrow(sig.assigns)
    num.loci<-as.numeric(gsub(file.patt,"\\1",filename))
    dir<-gsub("(.*)/.*.csv","\\1",filename)
    out<-data.frame(MarkerType=dir,NumLoci=as.numeric(num.loci),PropCorrect=as.numeric(prop.corr))
    rownames(out)<-filename
    return(out)
  },name.patt=name.patt,file.patt=file.patt))
  return(prop)
}
#####PLOT CERVUS INFO #####

stats.files<-c(list.files(path="parentage_haplotypes",pattern="\\d+.maternity.txt",full.names=TRUE),
               list.files(path="parentage_biallelic",pattern="\\d+.maternity.txt",full.names=TRUE))

stats<-cervus_analysis(stats.files,".*/[A-z]+(\\d+)_\\d+.maternity.txt",match_name = "Mother alone (all offspring):")
stats<-stats[!is.na(stats$NumLoci),]

r<-as.list(by(stats,stats$MarkerType,function(stat){
  rt<-tapply(as.numeric(stat[,"AssignmentRate"]),as.factor(stat[,"NumLoci"]),mean,na.rm=TRUE)
}))
png("CervusStats.png",height=5,width=10,res=300, units="in")
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(3,3,1,0.2))
plot_delta(stats,cols=c("slategray1","steelblue"),borders=c("slategray3","steelblue4"),"topright")
plot_assignmentRate(stats,cols=c("slategray3","steelblue4"),haps.name="parentage_haplotypes",
                    snps.name="parentage_biallelic",r=r,y.max = 30)
legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray3","steelblue4"),
	pch=c(19,15),lwd=2,lty=c(1,1),bty='n')
mtext("Number of Loci", 1,outer=T)
dev.off()

# How consistent are maternity assignments?
## Focus on biallelic markers with 1600 SNPs
mat.files<-list.files(path="parentage_biallelic/",pattern="dradPruned1600_\\d+_maternity.csv",full.names = TRUE)
cervus_mats<-function(files){
  out<-do.call(rbind,lapply(files,function(filename){
    dat<-read.csv(filename,row.names = NULL)
    colnames(dat)<-colnames(dat)[-1]
    off.mom<-dat[dat$Trio.confidence=="*",c("Offspring.ID","Candidate.mother.ID")]
    off.mom<-cbind(off.mom,filename)
    return(off.mom)
  }))
  return(out)
}
mats<-cervus_mats(mat.files)
off.mat<-tapply(mats$Candidate.mother.ID,mats$Offspring.ID,function(moms){
  length(unique(moms)) #the number of moms per offspring
})

length(unique(mats$Candidate.mother.ID))

##### PATERNITY #####
pat.files<-c(list.files(path="parentage_haplotypes",pattern="\\d+.paternity.txt",full.names=TRUE),
               list.files(path="parentage_biallelic",pattern="\\d+.paternity.txt",full.names=TRUE))

pats<-cervus_analysis(pat.files,".*/[A-z]+(\\d+)_\\d+.*paternity.txt",match_name = "Father alone (mother unknown):")
pats<-pats[!is.na(pats$NumLoci),]

pr<-as.list(by(pats,pats$MarkerType,function(stat){
  rt<-tapply(as.numeric(stat[,"AssignmentRate"]),as.factor(stat[,"NumLoci"]),mean,na.rm=TRUE)
}))


# evaluate whether they match
csv.files<-c(list.files(path="parentage_haplotypes",pattern="\\d+.paternity.csv",full.names=TRUE),
             list.files(path="parentage_biallelic",pattern="\\d+.paternity.csv",full.names=TRUE))

pat.corr<-rbind(cervus_patcorr(list.files(path="parentage_haplotypes",
                                    pattern="\\d+.paternity.csv",full.names=TRUE),
                         file.patt = ".*/[A-z]+(\\d+)_\\d+.*paternity.csv",
                         name.patt = "sample_\\w{3}([0-9]+)"),
                cervus_patcorr(list.files(path="parentage_biallelic",
                                    pattern="\\d+.paternity.csv",full.names=TRUE),
                         file.patt = ".*/[A-z]+(\\d+)_\\d+.*paternity.csv",
                         name.patt = "\\w+_(\\d+)"))


pat.assign<-as.list(by(pat.corr,pat.corr$MarkerType,function(stat){
  rt<-tapply(as.numeric(stat[,"PropCorrect"]),as.factor(stat[,"NumLoci"]),mean,na.rm=TRUE)
}))

png("Cervus_paternity.png",height=5,width=10,res=300, units="in")
par(mfrow=c(1,2),oma=c(1,1,1,1),mar=c(3,3,1,0.2))
#plot_delta(pats,cols=c("slategray1","steelblue"),borders=c("slategray3","steelblue4"),"topright")
plot_assignmentRate(pats,cols=c("slategray3","steelblue4"),haps.name="parentage_haplotypes",
                    snps.name="parentage_biallelic",r=pr)
#legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
#       pch=c(19,22),lwd=2,lty=c(1,1),bty='n')
plot_assignmentRate(pat.corr,cols=c("slategray3","steelblue4"),
                    haps.name="parentage_haplotypes",
                    snps.name="parentage_biallelic",r=pat.assign,
                    col2plot="PropCorrect",ylab="Proportion Correct Assignments",y.max = 1)
legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
       pch=c(19,15),lwd=2,lty=c(1,1),bty='n')
mtext("Number of Loci", 1,outer=T)
dev.off()

# what about only those with known dads?
KnownPat<-read.delim("colony/knownPat.txt",header=FALSE)
biKnownOff<-gsub("(\\w{2})\\w(\\d+.*)","\\1_\\2",KnownPat$V2)
bi.corr<-cervus_patcorr(list.files(path="parentage_biallelic",
                                   pattern="\\d+.paternity.csv",full.names=TRUE),
                        file.patt = ".*/[A-z]+(\\d+)_\\d+.*paternity.csv",
                        name.patt = "\\w+_(\\d+)",off.keep = biKnownOff)  
bicorr.means<-tapply(bi.corr$PropCorrect,bi.corr$NumLoci,mean)
  
##### SIMULATION OUTPUT #####
sim.files<-list.files(path="parentsim",pattern="_par.txt",full.names=TRUE)

sims<-cervus_analysis(sim.files,".*L(\\d+).*")
sims$MarkerType<-as.numeric(gsub(".*S(\\d+).*","\\1",sims$MarkerType))
sims$NumFemales<-as.numeric(gsub(".*F(\\d+).*","\\1",rownames(sims)))

r<-as.list(by(sims,sims$MarkerType,function(stat){
  rts<-by(stat,stat[,"NumFemales"],function(st){
    rt<-tapply(as.numeric(st[,"AssignmentRate"]),as.factor(st[,"NumLoci"]),mean,na.rm=TRUE)
  })
  
}))

par(mfrow=c(3,2),oma=c(1,1,1,1),mar=c(3,3,2,0.2))
#50 females
plot_delta(sims[sims$NumFemales==50 & !is.na(sims$Delta),],cols=c("slategray1","steelblue"),
           borders=c("slategray3","steelblue4"),"bottomleft")
plot_assignmentRate(sims[sims$NumFemales==50& !is.na(sims$Delta),],cols=c("slategray3","steelblue4"),
                    haps.name="4",snps.name="1",r=c(r[[2]][1],r[[1]][1]))
legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
       pch=c(19,22),lwd=2,lty=c(1,1),bty='n')
mtext("50 females",3,outer=TRUE,line=-0.5)
#100 females
plot_delta(sims[sims$NumFemales==100 & !is.na(sims$Delta),],cols=c("slategray1","steelblue"),
           borders=c("slategray3","steelblue4"),"bottomleft")
plot_assignmentRate(sims[sims$NumFemales==100& !is.na(sims$Delta),],cols=c("slategray3","steelblue4"),
                    haps.name="4",snps.name="1",r=c(r[[2]][2],r[[1]][2]))
legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
       pch=c(19,22),lwd=2,lty=c(1,1),bty='n')
mtext("100 females",3,outer=TRUE,line=-14.5)
#500 females
plot_delta(sims[sims$NumFemales==500 & !is.na(sims$Delta),],cols=c("slategray1","steelblue"),
           borders=c("slategray3","steelblue4"),"topleft")
plot_assignmentRate(sims[sims$NumFemales==500& !is.na(sims$Delta),],cols=c("slategray3","steelblue4"),
                    haps.name="4",snps.name="1",r=c(r[[2]][3],r[[1]][3]))
legend("bottomright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
       pch=c(19,22),lwd=2,lty=c(1,1),bty='n')
mtext("Number of Loci", 1,outer=T)
mtext("500 females",3,outer=TRUE,line=-29)

#####CHECK ALLELE FREQS
setwd("parentage_biallelic")
af.files<-list.files(pattern="\\d+.afs.txt")#or _allelefreqs
af.dat<-list()
for(i in 1:length(af.files)){
	nloci<-as.numeric(gsub("dradPruned(\\d+)_\\d+.afs.txt","\\1",af.files[i]))
	af<-read.table(af.files[i],skip=10,header=T,nrows=nloci)
	af.dat[[i]]<-as.data.frame(af)
}
afsum<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ cbind(mean(x$k),mean(x$HObs)) })))
colnames(afsum)<-c("MeanK","MeanObsHet")
rownames(afsum)<-af.files
afsum$NumLoci<-as.numeric(
	gsub("dradPruned(\\d+)_\\d+.afs.txt","\\1",rownames(afsum)))
afsum$setID<-gsub("dradPruned\\d+_(\\d+).afs.txt","\\1",rownames(afsum))
plot(afsum$NumLoci, afsum$MeanObsHet, xaxt='n',las=1,
	ylab="",xlab="",type='n',ylim=c(0,0.26))
text(labels=afsum$setID,x=afsum$NumLoci,y=afsum$MeanObsHet)
text(labels=stats$setID,x=stats$NumLoci,y=as.numeric(stats$AssignmentRate)/100,col="blue")
axis(1, at=c(50,100,150,200,300,400,800,1600),las=2)

nalleles<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ summary(x$k) })))
rownames(nalleles)<-af.files
nalleles$NumLoci<-as.numeric(
	gsub("dradPruned(\\d+)_\\d+.afs.txt","\\1",rownames(nalleles)))
nalleles$setID<-gsub("dradPruned\\d+_(\\d+).afs.txt","\\1",rownames(nalleles))
tapply(nalleles$Mean,nalleles$NumLoci,summary)

obshet<-data.frame(do.call("rbind",
	lapply(af.dat,function(x){ summary(x$HObs) })))
rownames(obshet)<-af.files
obshet$NumLoci<-as.numeric(
	gsub("dradPruned(\\d+)_\\d+.afs.txt","\\1",rownames(obshet)))
obshet$setID<-gsub("dradPruned\\d+_(\\d+)_afs.txt","\\1",rownames(obshet))
tapply(obshet$Mean,obshet$NumLoci,summary)

#####FULL DATASET
setwd("../")
full.af<-read.table("parentage/dradPrunedHaps_afs.txt",
	skip=10,header=T,nrows=124)#1642 for haplotypes

full.dat<-read.table("parentage/dradPrunedHaps.txt",header=T)
ids<-as.character(full.dat$ID)
pairs<-data.frame(id1=character(),id2=character())
for(i in 1:(length(ids)-1)){
	ida<-rep(ids[i],(length(ids)-i))
	idb<-ids[(i+1):length(ids)]
	pairs<-rbind(pairs,cbind(ida,idb))
}
write.table(pairs,"relatedness/pairwise.combinations.txt",col.names=F,row.names=F,
	quote=F,sep='\t')

gen.keep<-read.table("parentage_biallelic/PolymorphicIn90PercIndsHWE.txt",header=T)

obs.het<-NULL
sequence<-seq(2,ncol(gen.keep),2)
for(i in 1:length(sequence)){
	locus<-gen.keep[,c(sequence[i],sequence[i]+1)]
	loc<-locus[locus[,1] != "0",]
	hets<-nrow(loc[as.character(loc[,1]) != as.character(loc[,2]),])
	obs.het<-c(obs.het, (hets/nrow(loc)))
}
##############################################################################
####OTHER ANALYSES  parentsim_L400S1F500_2.crv
##############################################################################
setwd("../parentage")
hap.maternity.files<-list.files(pattern="\\d+_maternity.csv")
hap.maternity.dat<-data.frame()
for(i in 1:length(hap.maternity.files)){
	dat<-read.csv(hap.maternity.files[i],skip=1,header=F)
	colnames(dat)<-scan(hap.maternity.files[i],nlines=1,sep=",",
		what="character")
	sig<-dat[dat$"Trio confidence"=="*",]
	sig<-sig[,c("Offspring ID", "Candidate mother ID")]
	sig$NumLoci<-gsub("[A-z]+(\\d+)_\\d+.maternity.csv","\\1",
		hap.maternity.files[i])
	hap.maternity.dat<-rbind(hap.maternity.dat, sig)
}
hap.mat.split<-split(hap.maternity.dat, 
	factor(hap.maternity.dat$"Candidate mother ID"))
hap.summ.dat<- do.call("rbind", lapply(hap.mat.split,function(x){ 
	sum<-table(factor(x$"Offspring ID")) 
	min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ min(as.numeric(x$NumLoci)) }))
	max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ max(as.numeric(x$NumLoci)) }))
	y<-data.frame(MomID=levels(factor(x[,"Candidate mother ID"])),OffID=names(sum),
		NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc)
	rownames(y)<-NULL
	return(y)
}))
rownames(hap.summ.dat)<-NULL
write.csv(hap.summ.dat, "Cervus_Summary_Hap.csv")

##split by locus instead of females
hap.mat.loc<-split(hap.maternity.dat, factor(hap.maternity.dat$NumLoci))
hap.mat.loc.sum<-lapply(hap.mat.loc,function(ldf){ 
	df<-split(ldf,factor(ldf$"Candidate mother ID"))
	do.call("rbind", lapply(df, function(x){
		sum<-table(factor(x$"Offspring ID")) 
		min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ min(as.numeric(x$NumLoci)) }))
		max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ max(as.numeric(x$NumLoci)) }))
		y<-data.frame(MomID=levels(factor(x[,2])),OffID=names(sum),
			NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc)
		rownames(y)<-NULL
		return(y)
	}))
})
setwd("../parentage_biallelic")
snp.maternity.files<-list.files(pattern="\\d+.maternity.csv")
snp.maternity.dat<-data.frame()
for(i in 1:length(snp.maternity.files)){
	dat<-read.csv(snp.maternity.files[i],skip=1,header=F)
	colnames(dat)<-scan(snp.maternity.files[i],nlines=1,sep=",",
		what="character")
	sig<-dat[dat$"Trio confidence"=="*",]
	sig<-sig[,c("Offspring ID", "Candidate mother ID")]
	sig$NumLoci<-gsub("[A-z]+(\\d+)_\\d+.maternity.csv","\\1",snp.maternity.files[i])
	snp.maternity.dat<-rbind(snp.maternity.dat, sig)
}
snp.mat.split<-split(snp.maternity.dat, 
	factor(snp.maternity.dat$"Candidate mother ID"))
snp.summ.dat<- do.call("rbind", lapply(snp.mat.split,function(x){ 
	sum<-table(factor(x$"Offspring ID")) 
	min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ min(as.numeric(x$NumLoci)) }))
	max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
		function(x){ max(as.numeric(x$NumLoci)) }))
	y<-data.frame(MomID=levels(factor(x[,2])),OffID=names(sum),
		NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc,stringsAsFactors = FALSE)
	rownames(y)<-NULL
	return(y)
}))
rownames(snp.summ.dat)<-NULL
write.csv(snp.summ.dat, "Cervus_Summary_SNP.csv")

##split by locus instead of females
snp.mat.loc<-split(snp.maternity.dat, factor(snp.maternity.dat$NumLoci))

snp.mat.loc.sum<-lapply(snp.mat.loc,function(ldf){ 
	df<-split(ldf,factor(ldf$"Candidate mother ID"))
	do.call("rbind", lapply(df, function(x){
		sum<-table(factor(x$"Offspring ID")) 
		min.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ min(as.numeric(x$NumLoci)) }))
		max.loc<-unlist(lapply(split(x,factor(x$"Offspring ID")),
			function(x){ max(as.numeric(x$NumLoci)) }))
		y<-data.frame(MomID=levels(factor(x[,2])),OffID=names(sum),
			NumAssignments=sum,MinLoci=min.loc, MaxLoci=max.loc)
		rownames(y)<-NULL
		return(y)
	}))
})


###I'm not really using this.
png("CervusMaternitySummary.png",height=5,width=7,units="in",res=300)
par(mfrow=c(2,3), mar=c(2,2,2,2),oma=c(2,2,2,2))
hist(snp.mat.loc.sum$`100`$NumAssignments.Freq, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),ylim=c(0,130),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(hap.mat.loc.sum$`100`$NumAssignments.Freq, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),ylim=c(0,130),
	col=alpha("slategray4",0.5),border="slategray4",add=T)
legend("top",c("100 Loci"),bty='n')
mtext("Frequency", 2, outer=F,line=2,cex=0.75)
axis(1,pos=0)
axis(2, pos=0,las=1)
hist(snp.mat.loc.sum$`200`$NumAssignments, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(hap.mat.loc.sum$`200`$NumAssignments, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),
	col=alpha("slategray4",0.5),border="slategray4",add=T)
axis(1,pos=0)
axis(2, pos=0,las=1)
legend("top",c("200 Loci"),bty='n')
mtext("Number of Assignments", 1, outer=F,line=2,cex=0.75)
hist(snp.mat.loc.sum$`1600`$NumAssignments, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(hap.mat.loc.sum$`1600`$NumAssignments, breaks=seq(0,11,0.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",xlim=c(0,11),
	col=alpha("slategray4",0.5),border="slategray4",add=T)
axis(1,pos=0)
axis(2, pos=0,las=1)
legend("top",c("1600 Loci"),bty='n')

hist(summary(snp.mat.loc.sum$`100`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(summary(hap.mat.loc.sum$`100`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),add=T,
	col=alpha("slategray4",0.5),border="slategray4")
mtext("Frequency",2,outer=F,line=2,cex=0.75)
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
hist(summary(snp.mat.loc.sum$`200`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(summary(hap.mat.loc.sum$`200`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),add=T,
	col=alpha("slategray4",0.5),border="slategray4")
mtext("Inferred Mating Success",1,outer=F,line=2,cex=0.75)
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
hist(summary(snp.mat.loc.sum$`1600`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),
	col=alpha("steelblue4",0.5),border="steelblue4")
hist(summary(hap.mat.loc.sum$`1600`$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),
	xaxt='n',yaxt='n',main="",xlab="",ylab="",ylim=c(0,30),add=T,
	col=alpha("slategray4",0.5),border="slategray4")
axis(1,pos=0,at=c(1,2,3,4))
axis(2,pos=0.5, las=1)
legend("topright",c("Haplotypes","SNPs"),col=c("slategray4","steelblue4"),
	bty='n',pch=22,pt.bg=c(alpha("slategray4",0.5),alpha("steelblue4",0.5)))
dev.off()




###Old graph, don't use
#png("CervusMaternitySummary.png",height=7,width=7,units="in",res=300)
#par(mfrow=c(2,2),oma=c(2,2,2,2),mar=c(3,3,1,1))
#boxplot(summ.dat$NumAssignments~summ.dat$MaxLoci,las=1,ylab="",xlab="")
#mtext("Number of Assignments",2,outer=F,line=2,cex=0.75)
#mtext("Maximum Number of Loci",1,outer=F, line=2,cex=0.75)

#boxplot(summ.dat$NumAssignments~summ.dat$MinLoci,las=1,ylab="",xlab="")
#mtext("Minimum Number of Loci",1,outer=F, line=2,cex=0.75)

#hist(summary(summ.dat$MomID),breaks=c(0.5,1.5,2.5,3.5,4.5),xaxt='n',yaxt='n',
#	main="",xlab="",ylab="")
#mtext("Number of Offspring Assigned",1,outer=F,line=2,cex=0.75)
#mtext("Frequency",2,outer=F,line=2,cex=0.75)
#axis(1,pos=0,at=c(1,2,3,4))
#axis(2,pos=0.5, las=1)

#hist(summ.dat$NumAssignments, breaks=20,xaxt='n',yaxt='n',
#	main="",xlab="",ylab="")
#mtext("Number of Assignments", 1, outer=F,line=2,cex=0.75)
#axis(1,pos=0)
#axis(2, pos=0,las=1)
#dev.off()
