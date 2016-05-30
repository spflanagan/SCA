#Date: 1 February 2016
#Author: Sarah P. Flanagan
#Purpose: Choose a subsample of loci from a CERVUS genotypes file
#in which each locus is found in 98% of the individuals.
#then write that to file.

rm(list=ls())
#############################FUNCTIONS#####################################
prune.loci<-function(df, prop){
	num<-apply(df,2,function(x){ length(which(x=='0')) })
	keep<-num[num <= round(prop*nrow(df))]
	new.df<-df[,colnames(df) %in% names(keep)]
	return(new.df)
}

hwe.test<-function(df){#must be a df with two columns
	locus.name<-substr(names(df)[1],1,nchar(names(df)[1])-1)
	locus<-data.frame(A=df[df[,1]!="0",1],B=df[df[,2]!="0",2])
	gt.sort<-apply(locus,1,sort)
	locus$AB<-paste(gt.sort[1,],gt.sort[2,],sep="/")
	af<-table(c(as.character(locus$A),as.character(locus$B)))/
		sum(table(c(as.character(locus$A),as.character(locus$B))))
	gf<-table(locus$AB)
	result<-data.frame(genotype=character(),exp=numeric(),obs=numeric(),
		stringsAsFactors=F)
	for(i in 1:length(af)){ #calculate expected values
		for(j in 1:length(af)){
			names.sort<-sort(names(af)[c(i,j)])
			pq.name<-paste(names.sort[1],names.sort[2],sep="/")
			if(i == j){ 
				pq.exp<-af[i]*af[j]
			} else {	
				pq.exp<-2*af[i]*af[j]
			}
			if(length(gf[names(gf) %in% pq.name])>0){
				pq.obs<-gf[names(gf) %in% pq.name] 
				result[nrow(result)+1,]<-
					c(genotype=as.character(pq.name),exp=pq.exp,obs=pq.obs)
			} else {
				result[nrow(result)+1,]<-
					c(genotype=as.character(pq.name),exp=pq.exp,obs=rep(0,1))
			}
		}#end j
	}#end i
	#result has some duplicates because the order of alleles is not guaranteed.
	#test to see if expected = observed
	result<-result[!duplicated(result$genotype),]
	result$exp<-round(as.numeric(result$exp)*sum(as.numeric(result$obs)))
		result<-result[result$exp!="0",]
	result$chi<-((as.numeric(result$obs)-as.numeric(result$exp))^2)/
		as.numeric(result$exp)
	chi.result<-1-pchisq(sum(result$chi),length(af)-1)
	
	return(cbind(locus.name,chi.result))
}#end hwe.test function

#############################THE ACTUAL WORK################################
setwd("E://ubuntushare//SCA//results//parentage")
#large file
genotypes<-read.delim("genotypes.txt")
#genotypes<-read.delim("genotypes99_10loci.txt")

#prune to remove non-polymorphic loci
consensus.num<-apply(genotypes,2,function(x){ length(which(x=="consensus")) })
all.poly<-consensus.num[consensus.num==0]
poly<-genotypes[,colnames(genotypes) %in% names(all.poly)]

#remove any found in less than 90% of individuals
gen90<-prune.loci(poly,0.1)
write.table(gen90,"PolymorphicIn90PercInds.txt", quote=F,sep="\t",row.names=F)
#remove those loci not in hardy weinberg equilibrium
hwe<-data.frame()
for(x in seq(2,ncol(gen90),2)){
	hwe<-rbind(hwe,hwe.test(gen90[,c(x,(x+1))]))
}
keep.loci<-hwe[as.numeric(as.character(hwe$chi.result))>0.05,]
keep.names<-c("ID",unlist(lapply(as.list(keep.loci$locus.name),function(x){
	names<-c(paste(x,"A",sep=""),paste(x,"B",sep=""))
	return(names) })))
gen.keep<-gen90[,colnames(gen90) %in% keep.names]
write.table(gen.keep,"PolymorphicIn90PercIndsHWE.txt", quote=F,sep="\t",row.names=F)

#prune for higher coverage
pruned<-prune.loci(gen.keep, 0.01)
write.table(pruned,"PolymorphicIn99PercIndsHWE.txt", quote=F,sep="\t",row.names=F)
pruned<-read.table("PolymorphicIn99PercIndsHWE.txt")
#sample
nloci<-c(50,100,150,300,200,400,800,1600)
for(i in 1:10){ #do ten replicates of each set
	for(j in 1:length(nloci)){
		cols<-sample(seq(2,ncol(pruned),2),nloci[j],replace=F)
		cols<-c(cols,cols+1)
		write.table(pruned[,c(1,sort(cols))],
			paste("gen",nloci[j],"_",i,".txt",sep=""),
			quote=F,col.names=T,row.names=F,sep="\t")
	}
}

#######################PREVIOUS PRUNING
#multiple files
#genotypes.files<-list.files(pattern="genotypes[0-9]+.txt")
#all.pruned<-NULL
#for(i in 1:length(genotypes.files)){
#	gen<-read.delim(genotypes.files[i])
#	name<-sub("(genotypes[0-9]+).txt", "\\1", genotypes.files[i])
#	gen.new<-prune.loci(gen, 0.02)	
#	if(i == 1){
#		all.pruned<-data.frame(ID=gen[,1]) }
#	if(length(gen.keep) > 1){
#		all.pruned<-cbind(all.pruned,gen.new[,2:ncol(gen.new)]) 
#		print(paste(name," had ",ncol(gen.new)-1," present in 98% of individuals"))

		#write.table(gen.new,paste(name,".pruned.txt",sep=""), col.names=T,
			#row.names=F, sep='\t', quote=F) 
#	} else {
#		print(paste(name," had no loci present in 98% of individuals"))
#	}
#}

#write.table(gen90[,1],"CERVUS_indivduals.txt",row.names=F,sep='\t',quote=F,col.names=F)
#poly99<-prune.loci(poly, 0.01)

#write to file
#write.table(poly99,"genotypes_in99.txt", col.names=T,
#	row.names=F, sep='\t', quote=F)

#convert to relatedness format
#relatedness<-poly[,!(colnames(poly) %in% colnames(poly99))] #removes ID
#relatedness<-cbind(poly99, relatedness[,1:(2000-(ncol(poly99)-1))])
#write.table(relatedness,"relatedness1000.txt", col.names=T,
#	row.names=F, sep='\t', quote=F)