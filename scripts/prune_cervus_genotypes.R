#Date: 1 February 2016
#Author: Sarah P. Flanagan
#Purpose: Choose a subsample of loci from a CERVUS genotypes file
#in which each locus is found in 98% of the individuals.
#then write that to file.

rm(list=ls())
prune.loci<-function(df, prop){
	num<-apply(df,2,function(x){ length(which(x=='0')) })
	keep<-num[num <= round(prop*nrow(df))]
	new.df<-df[,colnames(df) %in% names(keep)]
	return(new.df)
}

setwd("E://ubuntushare//SCA//results//parentage")
#multiple files
genotypes.files<-list.files(pattern="genotypes[0-9]+.txt")
all.pruned<-NULL
for(i in 1:length(genotypes.files)){
	gen<-read.delim(genotypes.files[i])
	name<-sub("(genotypes[0-9]+).txt", "\\1", genotypes.files[i])
	gen.new<-prune.loci(gen, 0.02)	
	if(i == 1){
		all.pruned<-data.frame(ID=gen[,1]) }
	if(length(gen.keep) > 1){
		all.pruned<-cbind(all.pruned,gen.new[,2:ncol(gen.new)]) 
		print(paste(name," had ",ncol(gen.new)-1," present in 98% of individuals"))

		#write.table(gen.new,paste(name,".pruned.txt",sep=""), col.names=T,
			#row.names=F, sep='\t', quote=F) 
	} else {
		print(paste(name," had no loci present in 98% of individuals"))
	}
}

#large file
genotypes<-read.delim("genotypes.txt")
gen90<-prune.loci(genotypes,0.1)	
locus<-data.frame(A=gen90[gen90[,2]!="0",2],B=gen90[gen90[,3]!="0",3],
	AB=paste(gen90[gen90[,2]!="0",2],gen90[gen90[,3]!="0",3],sep="/"))
af<-table(c(as.character(locus$A),as.character(locus$B)))/
	sum(table(c(as.character(locus$A),as.character(locus$B)))
gf<-table(locus$AB)/sum(table(locus$AB))
result<-data.frame(genotype=character(),exp=numeric(),obs=numeric(),
	stringsAsFactors=F)
for(i in 1:length(af)){ #calculate expected values
	for(j in 1:length(af)){
		pq.name<-paste(names(af)[i],names(af)[j],sep="/")
		pq.exp<-2*af[i]*af[j]
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

hwt<-chisq.test(result$obs,result$exp)
	

###remove consensus
consensus.num<-apply(all.pruned,2,function(x){ length(which(x=="consensus")) })
all.poly<-consensus.num[consensus.num==0]
poly<-all.pruned[,colnames(all.pruned) %in% names(all.poly)]

pruned99<-prune.loci(all.pruned, 0.01)
poly99<-prune.loci(poly, 0.01)
write.table(poly99,"genotypes_in99.txt", col.names=T,
	row.names=F, sep='\t', quote=F)

relatedness<-poly[,!(colnames(poly) %in% colnames(poly99))] #removes ID
relatedness<-cbind(poly99, relatedness[,1:(2000-(ncol(poly99)-1))])
write.table(relatedness,"relatedness1000.txt", col.names=T,
	row.names=F, sep='\t', quote=F)

write.table(poly,"genotypes.poly.txt", col.names=T,
	row.names=F, sep='\t', quote=F) #all polymorphic loci
