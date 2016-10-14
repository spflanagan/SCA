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

fst.plot<-function(fst.dat,ci.dat, sig.col=c("red","yellow"),
	fst.name="Fst", chrom.name="Chrom", bp.name="BP",axis.size=0.5,
	scaffold.order=NULL,groups=NULL,print.names=FALSE,y.lim=NULL){
	if(!is.null(scaffold.order)){
		scaff.ord<-scaffold.order$component_id
		lgs<-scaffold.order$object
	} else{
		scaff.ord<-levels(factor(fst.dat[,chrom.name]))
		lgs<-scaff.ord
	}
	if(!is.null(groups)){
		lgs<-groups
		scaff.ord<-groups
	}
  fst.dat[,fst.name]<-as.numeric(as.character(fst.dat[,fst.name]))
  
	all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	xlist<-NULL
	xs<-NULL
	for(i in 1:length(scaff.ord)){
		all.scaff[[scaff.ord[i]]]<-
			all.scaff[[scaff.ord[i]]][order(all.scaff[[scaff.ord[i]]][,bp.name]),]	
		all.scaff[[scaff.ord[i]]][,bp.name]<-
			seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1)
		xs<-c(xs, seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1))
		new.max<-max(xs)
		#scaffold.order[i,"new_start"]<-last.max
		#scaffold.order[i,"new_end"]<-new.max
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		rownames(rect.xs)[i]<-scaff.ord[i]
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	#change BP to plot
	x.max<-max(xs)
	x.min<-min(xs)
	if(is.null(y.lim)){
  	y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
  	y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
  	if(min(fst.dat[,fst.name]) < 0) {
  		y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
  	} else {
  		y.min<-0
  	}
  
  	y.lim<-c(y.min,y.max)
	}
	displacement<-y.lim[1]-((y.lim[2]-y.lim[1])/30)
	plot(c(x.min,x.max),y.lim,xlim=c(x.min,x.max), 
		ylim=y.lim, 
		bty="n",type="n",	axes=F, xlab="", ylab="")
	for(i in 1:nrow(rect.xs)){
		if(i%%2 == 0) {
			rect.color<-"white"
		} else {
			rect.color<-"gray75"
		}
		rect(rect.xs[i,1],y.lim[1],rect.xs[i,2],y.lim[2], 
			col=rect.color, border=NA)
		if(print.names==T){
			text(x=mean(all.scaff[[scaff.ord[i]]][
				all.scaff[[scaff.ord[i]]]$Chrom==rownames(rect.xs)[i],
				bp.name]),
				y=displacement,labels=rownames(rect.xs)[i],
				adj=1,xpd=T,srt=45)
		}
	}
	for(i in 1:length(scaff.ord)){
		points(all.scaff[[scaff.ord[i]]][,bp.name], 
			all.scaff[[scaff.ord[i]]][,fst.name], 
			pch=19, cex=0.5,col="grey7",
			xlim=c(x.min,x.max),ylim=y.lim)
		#plot.genome.wide(all.scaff[[i]][,bp.name], 
		#	all.scaff[[i]][,fst.name],plot.rect=FALSE,
		#	y.max,x.max, y.min=y.min,x.min=x.min, 
		#	pt.col="grey7",#rect.xs[i,],rect.color,
		#	plot.new=TRUE, plot.axis=FALSE,  pt.cex=0.5)
		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] >= ci.dat[1],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col=sig.col[1], pch=19, cex=0.5)
		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] <= ci.dat[2],]
		points(temp.sig[,bp.name], temp.sig[,fst.name], 
			col=sig.col[2], pch=19, cex=0.5)
	}
	if(axis.size>0){
		axis(2, at = seq(round(y.lim[1],2),round(y.lim[2],2),
			round((y.lim[2]-y.lim[1])/2, digits=2)),
			ylim =y.lim, pos=0,
			labels=seq(round(y.lim[1],2),round(y.lim[2],2),
				round((y.lim[2]-y.lim[1])/2, digits=2)),
			las=1,tck = -0.01, xlab="", ylab="", cex.axis=axis.size)
	}
	xes<-do.call("rbind",all.scaff)
	return(xes)
}

