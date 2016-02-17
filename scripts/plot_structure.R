

plot.structure<-function(structure.out, k, pop.order, 
	filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,
	plot.new=TRUE){
	str.split<-split(structure.out,structure.out[,1])
	bar.colors<-rainbow(k,s=0.5)
	if(make.file==TRUE){
		jpeg(filename,width=7, height=1.25, units="in", res=300)
		par(mfrow=c(1,length(str.split)))
	} 
	#par(mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
	for(i in 1:length(str.split)){
		pop.index<-pop.order[i]
		barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
			beside=FALSE, space=0,	border=NA, col=bar.colors,
			xlab="", ylab="", xaxt='n', yaxt='n')#, new=plot.new)
		mtext(pop.index, 1, line=0.5, cex=0.5, outer=F)
	}
	if(make.file==TRUE) {dev.off()}
}

setwd("E://ubuntushare//SCA//results//structure//sca//admixture//Results//")
structure.files<-list.files(pattern="*f_clusters.txt")
structure<-lapply(structure.files,read.table)
structure<-lapply(structure, function(x){
	x$V1<-sub('sample_([A-Z]{3})\\d','\\1',x$V1) 
	return(x)})
k<-sub('\\w+_\\w+_(\\d+)_f_clusters.txt','\\1',structure.files)
k[k %in% seq(1,10)]<-paste("1_",k[k %in% seq(1,10)],sep="")
k[k %in% seq(11,20)]<-paste("2_",k[k %in% seq(11,20)],sep="")
k[k %in% seq(21,30)]<-paste("3_",k[k %in% seq(21,30)],sep="")
k[k %in% seq(31,40)]<-paste("4_",k[k %in% seq(31,40)],sep="")
names(structure)<-k
pop.list<-levels(as.factor(sub('sample_([A-Z]{3})\\d','\\1', structure[[1]]$V1)))
	

#par(mfrow=c(4,length(pop.list)),mar=c(1,0,1,0),oma=c(1,0,1,0))
png("giant.structure.png",height=50,width=10,units="in",res=300)
par(mfrow=c(length(structure),length(pop.list)),mar=c(1,0,1,0),oma=c(1,0,1,0),cex=0.5)
plot.structure(structure[[1]],
		sub('(\\d+)_\\d+','\\1',names(structure[1])),
		pop.list, make.file=FALSE)#, plot.new=FALSE)
for(i in 2:length(structure)){
	plot.structure(structure[[i]],
		sub('(\\d+)_\\d+','\\1',names(structure[i])),
		pop.list, make.file=FALSE, plot.new=FALSE)

}
dev.off()

#just plot one example of each set.
png("sca.structure.png",height=10,width=10,units="in",res=300)
par(mfrow=c(3,length(pop.list)),mar=c(1,0,1,0),oma=c(1,0,1,0),cex=0.5)
plot.structure(structure[[11]],
		sub('(\\d+)_\\d+','\\1',names(structure[11])),
		pop.list, make.file=FALSE)#, plot.new=FALSE)
plot.structure(structure[[21]],
		sub('(\\d+)_\\d+','\\1',names(structure[21])),
		pop.list, make.file=FALSE)#, plot.new=FALSE)
plot.structure(structure[[31]],
		sub('(\\d+)_\\d+','\\1',names(structure[31])),
		pop.list, make.file=FALSE)#, plot.new=FALSE)
dev.off()



