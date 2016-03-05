#Author: Sarah P. Flanagan
#Date: 26 Januray 2016
#Purpose: To link sample*.snps.tsv file to catalog IDs and rewrite for gwsca

reformat.snp<-function(filename){
	matches<-read.table(paste("stacks/",filename,".matches.tsv",sep=""))
	colnames(matches)<-c("SQLID","BatchID","CatID","SampleID","LocusID",
		"Hap","StackDepth","LogLike")
	snps<-read.table(paste("stacks/",filename,".snps.tsv",sep=""))
	colnames(snps)<-c("SQLID","SampleID","LocusID","Column","Type","LR",
		"MajNuc","AltNuc")

	#remove those with likelihood ratios too low to call genotype
	snps<-snps[snps$Type!="U",]
	
	#merge to align different SNP IDs
	new.snp<-merge(snps[,c("LocusID","Column","MajNuc","AltNuc","Type")], 
		matches[,c("CatID","LocusID","StackDepth")], by="LocusID",all=F)
	new.snp$SNPID<-paste(new.snp$CatID,".",new.snp$Column,sep="")
	
	#compare to reference
	keep.snp<-merge(new.snp,inform.snps,by="SNPID")
	
	#prune to remove any with StackDepth < 5
	keep.snp<-keep.snp[keep.snp$StackDepth >= 5,]

	#convert majnuc and altnuc to based on Type
	keep.snp[which(keep.snp$Type=="O"),"AltNuc"]<-
		keep.snp[which(keep.snp$Type=="O"),"MajNuc"]

	out.snp<-keep.snp[,c("SNPID","MajNuc","AltNuc")]
	write.table(out.snp, paste("biallelic/",filename,".new.snps.txt",sep=""), 
		sep="\t", quote=F, col.names=T,row.names=F)
	return(ncol(out.snp))
}

setwd("E://ubuntushare//SCA//results//")

sumstats<-read.table("stacks/batch_1.sumstats.tsv")
#colnames(sumstats)<-c("batch","LocID","Chr","BP","Col","PopID","P","Q",
#	"N","PFreq","ObsHet","ObsHom","ExpHet","ExpHom","pi","spi","spiP",
#	"Fis","SFis","SFisP","Private")
#sumstats$SNPID<-paste(sumstats$LocID,".",sumstats$Col,sep="")


inform.snps<-sumstats[,2:5]
colnames(inform.snps)<-c("LocID","Chr","BP","Col")
inform.snps$SNPID<-paste(inform.snps$LocID,".",inform.snps$Col,sep="")
rm(sumstats)

files<-read.table("biallelic/samples_list.txt", stringsAsFactors=F)
reformat<-lapply(files$V1, reformat.snp)


#dads<-cbind(dads,gsub(".*(\\d{3}).*","\\1",dads))
#colnames(dads)<-c("id","num")
#off<-cbind(off,gsub(".*(\\d{3}).*","\\1",off))
#colnames(off)<-c("id","num")
#write.table(merge(dads,off,by="num"),"dad.kid.pairs.fullnames.txt",quote=F,
#	col.names=F,row.names=F,sep="\t")