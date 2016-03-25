
##TROUBLESHOOTING THE WEIRD PATTERNS
#CALCULATE FST HERE
fst.aj<-merge(adt.n,juv.n,by="Locus")
fst.aj$AvgHs<-(fst.aj$Hs.x+fst.aj$Hs.y)/2
fst.aj$pbar<-((fst.aj$Allele1Freq.x*fst.aj$N.x)+
	(fst.aj$Allele1Freq.y*fst.aj$N.y))/(fst.aj$N.x+fst.aj$N.y)
fst.aj$qbar<-((fst.aj$Allele2Freq.x*fst.aj$N.x)+
	(fst.aj$Allele2Freq.y*fst.aj$N.y))/(fst.aj$N.x+fst.aj$N.y)
fst.aj$ht<-1-((fst.aj$pbar*fst.aj$pbar)+(fst.aj$qbar*fst.aj$qbar))
fst.aj$fst<-(fst.aj$ht-fst.aj$AvgHs)/fst.aj$ht

##TROUBLESHOOTING MALE-FEMALE COMPARISON
 weirdos<-fm.prune[fm.prune$MAL.FEM > 0.03,]
weirdos<-weirdos[,c("Locus","Chrom","Pos","MAL.FEM")]
weirdsum<-gw.sum[gw.sum$Locus %in% weirdos$Locus,]
weirdsum<-weirdsum[weirdsum$Pop == "MAL" | weirdsum$Pop == "FEM",]
weirdos$hobs<-1-weirdos$Ho
regsum<-gw.sum[gw.sum$Locus %in% fm.prune$Locus,] 
regsum<-regsum[regsum$Pop=="MAL" | regsum$Pop == "FEM",]
regsum<-regsum[!(regsum$Locus %in% weirdsum$Locus),]
regsum$hobs<-1-regsum$Ho

png("HsvHo.png",height=7,width=7,units="in",res=300)
par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(regsum[regsum$Pop=="FEM",c("hobs","Hs")])
mtext("Regular, Females",3,outer=F)
plot(regsum[regsum$Pop=="MAL",c("hobs","Hs")])
mtext("Regular, Males",3,outer=F)
plot(weirdos[weirdos$Pop=="FEM",c("hobs","Hs")])
mtext("Weird, Females",3,outer=F)
plot(weirdos[weirdos$Pop=="MAL",c("hobs","Hs")])
mtext("Weird, Males",3,outer=F)
mtext("Observed Heterozygosity",1,outer=T)
mtext("Expected Heterozygosity (Hs)",2,outer=T)
dev.off()

png("HsInWeirdMalFemLoci.png",height=7,width=7,units="in",res=300)
par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
hist(weirdsum[weirdsum$Pop=="MAL","Hs"], breaks=50, main="",ylab="",xlab="")
legend("topright",bty="n","Weird Loci, Males")
hist(regsum[regsum$Pop=="MAL","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty="n","Regular Loci, Males")
hist(weirdsum[weirdsum$Pop=="FEM","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty='n',"Weird Loci, Females")
hist(regsum[regsum$Pop=="FEM","Hs"], breaks=50, main="",ylab="",xlab="")
legend("top",bty='n',"Regular Loci, Females")
mtext("Hs",1,outer=T)
mtext("Frequency",2,outer=T)
dev.off()


png("HvNInWeirdMalFemLoci.png",height=7,width=7,units="in",res=300)
par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(weirdsum[weirdsum$Pop=="MAL",c("N","Hs")],ylab="",xlab="",pch=15)
legend("bottomleft",bty="n","Hs Males",text.col="red")
plot(weirdsum[weirdsum$Pop=="MAL",c("N","Ho")],ylab="",xlab="",pch=15,col="blue")
legend("bottomleft",bty="n","Ho, Males",text.col="red")
plot(weirdsum[weirdsum$Pop=="FEM",c("N","Hs")],ylab="",xlab="",pch=19)
legend("bottomright",bty='n',"Hs, Females",text.col="red")
plot(weirdsum[weirdsum$Pop=="FEM",c("N","Ho")],ylab="",xlab="",pch=19,col="blue")
legend("bottomleft",bty='n',"Ho, Females",text.col="red")
mtext("Heterozygosity",2,outer=T)
mtext("N",1,outer=T)
mtext("Weird Loci",3,outer=T)
dev.off()

 par(mfrow=c(2,2), oma=c(2,2,2,2),mar=c(2,2,2,2))
plot(regsum[regsum$Pop=="MAL",c("N","Hs")],ylab="",xlab="",pch=15)
legend("bottomleft",bty="n","Hs Males",text.col="red")
plot(regsum[regsum$Pop=="MAL",c("N","Ho")],ylab="",xlab="",pch=15,col="blue")
legend("bottomleft",bty="n","Ho, Males",text.col="red")
plot(regsum[regsum$Pop=="FEM",c("N","Hs")],ylab="",xlab="",pch=19)
legend("bottomright",bty='n',"Hs, Females",text.col="red")
plot(regsum[regsum$Pop=="FEM",c("N","Ho")],ylab="",xlab="",pch=19,col="blue")
legend("bottomleft",bty='n',"Ho, Females",text.col="red")
mtext("Heterozygosity",2,outer=T)
mtext("N",1,outer=T)
mtext("Weird Loci",3,outer=T)


fem.n<-sum.list$FEM[!is.na(sum.list$FEM$Hs),]
 plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#OK
fem.n<-fem.n[fem.n$Allele1Freq > 0.05 & fem.n$Allele1Freq < 0.95,]
plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#OK
fem.n<-fem.n[fem.n$N>100,]
plot(gw.fst[gw.fst$Locus %in% fem.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])#goodbye

removedfem<-fem.n[fem.n$N < 100,]
plot(gw.fst[gw.fst$Locus %in% removedfem$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])
rem.fem<-gw.fst[gw.fst$Locus %in% removedfem$Locus & gw.fst$MAL.FEM > 0.025 
	& gw.fst$MAL.FEM < 0.05,c("Locus","Chrom","Pos","MAL.FEM")]
rem.fem<-merge(rem.fem, removedfem,by="Locus")

rem.mal<-gw.sum[gw.sum$Locus %in% rem.fem$Locus & gw.sum$Pop=="MAL",]

mal.n<-sum.list$MAL[sum.list$MAL$N>160& !is.na(sum.list$MAL$Hs),]
 plot(gw.fst[gw.fst$Locus %in% mal.n$Locus & gw.fst$MAL.FEM > 0,"MAL.FEM"])

png("CoverageOfLocusGroups.png",height=7,width=14,units="in",res=300)
par(mfrow=c(2,4), oma=c(2,2,2,2),mar=c(2,2,2,2))
hist(regsum[regsum$Pop=="MAL","N"], breaks=50, main="",ylab="",xlab="")
mtext("Regular (kept) Loci, Males",3,outer=F)
hist(weirdsum[weirdsum$Pop=="MAL","N"], breaks=50, main="",ylab="",xlab="")
mtext("Weird Loci, Males",3,outer=F)
hist(gw.sum[gw.sum$Locus %in% removedfem$Locus & gw.sum$Pop=="MAL","N"], breaks=50, main="",ylab="",xlab="")
mtext("All Pruned-Out Loci, Males",3,outer=F)
hist(rem.mal$N, breaks=50, main="",ylab="",xlab="")
mtext("Pruned-Out Loci in Fst Gap, Males",3,outer=F)
hist(regsum[regsum$Pop=="FEM","N"], breaks=50, main="",ylab="",xlab="")
mtext("Regular (kept) Loci, Females",3,outer=F)
hist(weirdsum[weirdsum$Pop=="FEM","N"], breaks=50, main="",ylab="",xlab="")
mtext("Weird Loci, Females",3,outer=F)
hist(removedfem$N, breaks=50, main="",ylab="",xlab="")
mtext("All Pruned-Out Loci, Females",3,outer=F)
hist(rem.fem$N, breaks=50, main="",ylab="",xlab="")
mtext("Pruned-Out Loci in Fst Gap, Females",3,outer=F)
mtext("N",1,outer=T)
mtext("Frequency",2,outer=T)
dev.off()

gw.uber<-merge(gw.sum,gw.fst,by="Locus") 
png("NvFst.png",height=7,width=7,units="in",res=300)
par(mfrow=c(1,2))
plot(gw.uber[gw.uber$Pop=="MAL" & gw.uber$MAL.FEM > 0,c("N","MAL.FEM")],ylab="Fst")
mtext("Males",3,outer=F)
plot(gw.uber[gw.uber$Pop=="FEM" & gw.uber$MAL.FEM > 0,c("N","MAL.FEM")],ylab="Fst")
mtext("Females",3,outer=F)
dev.off()



tags<-read.table("../stacks/batch_1.catalog.tags.tsv",sep='\t')
colnames(tags)<-c("SqlID","SampleID","LocusID","Chr","BP","Strand","SeqType",
	"StackComponent","SeqID","Sequence","Deleveraged","Blacklist",
	"Lumberjackstack","loglike")
snps<-read.table("../stacks/batch_1.catalog.snps.tsv",sep='\t')
colnames(snps)<-c("SqlID","SampleID","LocusID","Col","Type","LikelihoodRatio",
	"Rank1","Rank2","Rank3","Rank4")
sumstats<-read.table("../stacks/batch_1.sumstats.tsv",sep='\t')
colnames(sumstats)<-c("BatchID","LocusID","Chr","BP","Col","PopID","P","Q",
	"NumInd","FreqP","ObsHet","ObsHom","ExpHet","ExpHom","Pi","SmoothPi",
	"SmoothPiP","Fis","SmoothFis","SmoothFisP","Private")
sumstats$Locus<-paste(sumstats$Chr,".",sumstats$BP,sep="")
tagsnps<-merge(tags,snps,"LocusID")
tagsnps$BPCol<-tagsnps$BP+tagsnps$Col+1
tagsnps$Locus<-paste(tagsnps$Chr,".",tagsnps$BPCol,sep="")
tagsnps<-merge(tagsnps,sumstats,"Locus")
weird.tags<-tagsnps[tagsnps$Locus %in% weirdos$Locus,]
weird.tags<-weird.tags[,c("LocusID","Chr.x","BP.x","Chr.y","BP.y",
	"Sequence","loglike","LikelihoodRatio","Col.y","PopID","Rank1",
	"Rank2","P","Q","NumInd","FreqP","ObsHet","ObsHom","ExpHet",
	"ExpHom","Pi","Fis","Private","Locus")]
weird.tags.out<-merge(weirdos,tagsnps,"Locus")
write.table(weird.tags.out,"WeirdFstsInfo.txt",sep='\t',quote=F,row.names=F)
#...what now???