#Author: Sarah P. Flanagan
#Last Updated: 27 April 2016
#Start Date: 27 April 2016
#Purpose: Conduct a linkage disequilibrium analysis

genome<-read.delim("E:/ubuntushare/scovelli_genome/SSC_genome.agp",	
	comment.char="#",header=F)
scaffs<-read.delim("E:/ubuntushare/scovelli_genome/SSC_scaffolds.agp",	
	comment.char="#",header=F)
colnames(genome)<-c("object","object_beg","object_end","part_number",
	"W","component_id","component_beg","component_end","orientation")
colnames(scaffs)<-c("object","object_beg","object_end","part_number",
	"W","component_id","component_beg","component_end","orientation")

gen<-genome[genome$W=="W",]
scaffold.order<-gen$component_id
fst.gen<-merge(gen,gw.fst,by.x="component_id",by.y="Chrom",all.y=T)

fst.gen<-gw.fst[match(scaffold.order, gw.fst$Chrom), ]

gff<-read.delim("E:/ubuntushare/scovelli_genome/annotated_scovelli_scaffolds_March2016.gff",
	header=F)
gff.chrom<-levels(gff$V1)
gw.chrom<-levels(gw.fst$Chrom)


vcf<-read.delim("../stacks/batch_1.vcf",comment.char="#",sep='\t')
header<-scan("../stacks/batch_1.vcf",what="character")[header.start:
	(header.start+ncol(vcf)-1)]
colnames(vcf)<-header1
if(length(strsplit(as.character(vcf[1,10]),":")[[1]])>1){
	new<-vcf[,1:3]
	for(i in 10:ncol(vcf)){
	new<-cbind(new,
		sapply(vcf[,i],function(x) {
			strsplit(as.character(x),":")[[1]][1]})
		)
	}
	colnames(new)<-colnames(vcf[,c(1:3,10:ncol(vcf))])
	vcf<-new
}
colnames(vcf)<-header[c(1:3,10:length(header))]
vcf.chrom<-split(vcf,vcf$`#CHROM`)

calc.ld<-function(vcf.list,row1,row2){
	loc1<-vcf.list[row1,4:ncol(vcf.list)]
	loc2<-vcf.list[row2,4:ncol(vcf.list)]
	joint<-rbind(loc1,loc2)
	joint<-joint[,joint[1,] != "./." & joint[2,]!="./."]
	joint.freqs<-data.frame("0"=c(0,0),"1"=c(0,0),row.names=c("0","1"))
	freqs1<-data.frame("0"=0,"1"=0)
	freqs2<-data.frame("0"=0,"1"=0)
	for(i in 1:ncol(joint)){
		mat1<-as.numeric(strsplit(as.character(joint[1,i]),"/")[[1]][1])+1
		pat1<-as.numeric(strsplit(as.character(joint[1,i]),"/")[[1]][2])+1
		mat2<-as.numeric(strsplit(as.character(joint[2,i]),"/")[[1]][1])+1
		pat2<-as.numeric(strsplit(as.character(joint[2,i]),"/")[[1]][2])+1
		freqs1[mat1]<-freqs1[mat1]+1
		freqs1[pat1]<-freqs1[pat1]+1
		freqs2[mat2]<-freqs2[mat2]+1
		freqs2[pat2]<-freqs2[pat2]+1
		joint.freqs[mat1,mat2]<-joint.freqs[mat1,mat2]+1
		joint.freqs[pat1,pat2]<-joint.freqs[pat1,pat2]+1
	}
	joint.freqs<-joint.freqs/(2*ncol(joint))
	freqs1<-freqs1/(2*ncol(joint))
	freqs2<-freqs2/(2*ncol(joint))
	
	d<-data.frame("0"=c(0,0),"1"=c(0,0),row.names=c("0","1"))
	dmax<-d
	for(i in 1:2){
		for(j in 1:2){	
			d[i,j]<-joint.freqs[i,j]-(freqs1[i]*freqs2[j])
			if(d[i,j]<0){	
				dm<-min((freqs1[i]*freqs2[j]),
					((1-freqs1[i])*(1-freqs2[j])))
			}else{
				dm<-min(((1-freqs1[i])*freqs2[j]),
					((freqs1[i])*(1-freqs2[j])))
			}
			dmax[i,j]<-dm
		}
	}
	dprime<-0
	for(i in 1:2){
		for(j in 1:2){
			if(freqs1[i] > 0 & freqs2[j] > 0){
				if(dmax[i,j]>0){
					dprime<-dprime+(freqs[i]*freqs[j]*
						abs(d[i,j])/dmax[i,j])
				} else { dprime<- -5 }
			}
		}
	}
	return(dprime)
}

ld<-list(rep(data.frame(),length(vcf.chrom)))
for(f in 1:length(vcf.chrom)){
	for(ff in 1:length(vcf.chrom[[f]])){
		for(fff in 1:length(vcf.chrom[[f]])){
			ld[[f]][ff,fff]<-calc.ld(vcf.chrom[[f]],ff,fff)
		}
	}

}

snpstats<-read.table("../sexlinked/snpstats1_out.txt",header=T)
snp.plots<-snpstats[,c("LG","pos","p_val.3")]
snp.plots$plotp<--log10(snp.plots$p_val.3)
snp.plots<-snp.plots[snp.plots$plotp != "Inf" &snp.plots$plotp != "-Inf",]
png("logP_G.png", width=10,height=7,units="in",res=300)
g<-fst.plot(snp.plots,ci.dat=c(0,0),sig.col=c("black","black"),
	fst.name="plotp",chrom.name="LG",bp.name="pos")
dev.off()

#ldc<-read.table("ld_info.txt",header=T)
ldc<-read.table("ld_matrix_lg1.txt",header=T,row.names=1)
ldc[ldc==-1]<-NA

#heatmap.2(as.matrix(ldc),dendrogram="none",density.info='none',tracecol="NA")


ldc.split<-split(ldc,ldc$LocIDA)
ld.mat<-matrix(nrow=length(levels(factor(ldc$LocIDA))),
	ncol=length(levels(factor(ldc$LocIDB))),
	dimnames=list(levels(factor(ldc$LocIDA)),levels(factor(ldc$LocIDB))))
for(i in 1:nrow(ld.mat)){
	for(j in 1:nrow(ld.mat)){
		ld.mat[i,j]<-ldc[ldc$LocIDA==rownames(ld.mat)[i] &
			ldc$LocIDB==rownames(ld.mat)[j],"D."]
	}
}
starts<-c(seq(1,4581,20))
ends<-c(seq(20,4580,50),4603)
pdf("LDheatmap.LG1.pdf",height=10,width=10)
par(mfrow=c(length(starts),length(starts)),oma=c(0,0,0,0),mar=c(0,0,0,0))
for(i in 1:length(starts)){
	for(j in 1:length(ends)){
		#heatmap.2(as.matrix(ldc[starts[i]:ends[i],starts[j]:ends[j]]),
		#	dendrogram="none",tracecol="NA",labCol="",labRow="",key=F,
		#	Colv=F,Rowv=F,lwid=c(0.5,4),lhei=c(0.5,4),new=F)
		image(as.matrix(ldc[starts[i]:ends[i],starts[j]:ends[j]]),
			xaxt='n',yaxt='n')
		print(paste(starts[i],":",ends[i],",",starts[j],":",ends[j],sep=""))
}}
dev.off()




