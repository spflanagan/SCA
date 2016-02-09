#Author: Sarah P. Flanagan
#Date: 4 February 2016
#Purpose: Plot het_v_depth results from Monnahan et al. (2015)'s SCA

path<-"E:/ubuntushare/SCA/results/monnahan"
setwd(path)
hd<-read.delim("het_v_depth.txt",header=T)
hd$prop.het<-hd$called.het/(hd$called.het+hd$called.homo)
if(min(hd[which(rowSums(hd[,2:3])==0),"read.depth"])<100){
	x.max<-min(hd[which(rowSums(hd[,2:3])==0),"read.depth"])
} else {
	x.max<-100 }

png("het_v_depth.png",height=7,width=7, units="in",res=300)
plot(hd$read.depth,hd$prop.het,pch=19, xlab="Read Depth",
	ylab="Proportion Called Heterozygote",las=1, xlim=c(0,x.max))
dev.off()
