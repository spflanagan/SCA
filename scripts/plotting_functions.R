#Author: Sarah P. Flanagan 
#Purpose: Plot genome-wide statistics.

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plot.genome.wide<-function(bp,var,y.max,x.max, rect.xs=NULL,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white",plot.rect=TRUE, 
	pt.cex=1, pt.col="black"){
	#********************************************
	#this function plots a variable without scaffold info. 
	#feed it the basepair (x) values and variable (y) values 
	#*********************************************
	if(plot.new==TRUE){ par(new=new) }
	plot(bp, var,xlab="",ylab="", new=plot.new,
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	if(plot.rect==TRUE){
		num.rect<-nrow(rect.xs)
		if(is.null(num.rect)) {
			rect(rect.xs[1],y.min,rect.xs[2],y.max, 
				col=rect.color, border=NA)
		} else {
			for(i in 1:nrow(rect.xs)){
				rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
					col=rect.color, border=NA)
			}
		}
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,col=pt.col,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
#PLOT ALL THE SCAFFOLDS
#***************************************************************************#

plot.fsts<-function(fst.dat,ci.dat, 
	fst.name="Fst", chrom.name="Chrom", bp.name="BP",axis.size=0.5){
	all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	for(i in 1:length(all.scaff)){
		new.max<-last.max+round(max(all.scaff[[i]][,bp.name]), -2)
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	#change BP to plot
	for(i in 1:length(all.scaff)){
		all.scaff[[i]][,bp.name]<-
			all.scaff[[i]][,bp.name]+addition.values[i]
	}
	x.max<-max(addition.values)
	x.min<-min(all.scaff[[1]][,bp.name])
	y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
	y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
	if(min(fst.dat[,fst.name]) < 0) {
		y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
	} else {
		y.min<-0
	}

	plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), 
		ylim=c(y.min, y.max), 
		bty="n",type="n",	axes=F, xlab="", ylab="")
	for(i in 1:nrow(rect.xs)){
		if(i%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray96"
		}
		rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
			col=rect.color, border=NA)
	}
	for(i in 1:length(all.scaff)){
		plot.genome.wide(all.scaff[[i]][,bp.name], 
			all.scaff[[i]][,fst.name],plot.rect=FALSE,
			y.max,x.max, rect.xs[i,],y.min=0,x.min=0, pt.col="grey53",
			plot.new=TRUE, plot.axis=FALSE, rect.color, pt.cex=0.5)
		temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] >= ci.dat[1],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col="red", pch=19, cex=0.5)
		temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] <= ci.dat[2],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col="yellow", pch=19, cex=0.5)
	}
	axis(2, at = seq(round(y.min,2),round(y.max,2),
			round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		labels=seq(round(y.min,2),round(y.max,2),
			round((y.max-y.min)/2, digits=2)),
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=axis.size)
}


