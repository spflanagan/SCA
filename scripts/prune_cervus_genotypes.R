#Date: 1 February 2016
#Author: Sarah P. Flanagan
#Purpose: Choose a subsample of loci from a CERVUS genotypes file
#in which each locus is found in 98% of the individuals.
#then write that to file.

rm(list=ls())
#############################FUNCTIONS#####################################
prune.loci<-function(df, prop,blanks='0'){
	num<-apply(df,2,function(x){ length(which(x==blanks)) })
	keep<-num[num <= round(prop*nrow(df))]
	new.df<-data.frame(df[,colnames(df) %in% names(keep)])
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


#sample
rep_cervus<-function(nloci,nreps,gens,out.prefix="gen"){#used pruned for haplotypes, gen.keep for SNPs
  for(i in 1:nreps){ #do nreps replicates of each set
    for(j in 1:length(nloci)){
      cols<-sample(seq(2,ncol(gens)-2,2),nloci[j],replace=F)
      cols<-c(cols,cols+1)
      write.table(gens[,c(1,sort(cols))],
                  paste(out.prefix,nloci[j],"_",i,".txt",sep=""),
                  quote=F,col.names=T,row.names=F,sep="\t")
    }
  }
}

#############################THE ACTUAL WORK################################
setwd("B://ubuntushare//SCA//results//parentage_biallelic")
#large file
genotypes<-read.delim("snp_genotypes_batch_1.txt")
#genotypes<-read.delim("snp_genotypes.txt")
#genotypes<-read.delim("genotypes99_10loci.txt")


#remove any found in less than 90% of individuals
gen90<-prune.loci(genotypes,0.1)
write.table(gen90,"PolymorphicIn90PercInds.txt", quote=F,sep="\t",row.names=F)
#remove those loci not in hardy weinberg equilibrium
hwe<-data.frame()
for(x in seq(2,(ncol(gen90)-1),2)){
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
#only 121 snps now with drad_miss3.vcf 
write.table(pruned,"PolymorphicIn99PercIndsHWE.txt", quote=F,sep="\t",row.names=F)
pruned<-read.table("PolymorphicIn99PercIndsHWE.txt")

#For CERVUS
s<-rep_cervus(nloci=2000,nreps=1,gens=gen.keep,out.prefix = "dradPruned")
c<-rep_cervus(nloci=c(50,100,150,300,200,400,800,1600),nreps=10,gens=gen.keep,out.prefix = "dradPruned")
write.table(gen.keep$ID[grep("FE",gen.keep$ID)],"candidate_moms.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
off<-gen.keep$ID[grep("OF",gen.keep$ID)]
dad<-gen.keep$ID[grep("PR",gen.keep$ID)]
offs<-data.frame(offspring=off,KnownParent=NA,stringsAsFactors = FALSE)
for(i in 1:length(off)){
  d<-which(dadn==offn[i])
  if(length(d)==0){ offs[i,2]<-"" 
  }else { offs[i,2]<-as.character(dad[d]) }
}
write.table(offs,"offspring.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

#############################HAPLOTYPES################################
setwd("B://ubuntushare//SCA//results//parentage")

#convert haplotypes to cervus format
haps<-read.delim("../stacks/batch_1.haplotypes.tsv")
haps.miss<-prune.haps(haps,0.3)
haps.keep<-do.call(cbind,apply(haps.miss,1,function(hap){
  ids<-as.character(hap["Catalog.ID"])
  idsb<-paste(ids,"b",sep="")
  idsa<-paste(ids,"a",sep="")
  cnt<-hap["Cnt"]
  hap<-unlist(lapply(hap[3:length(hap)],as.character))
  hap[hap=="-"]<-"0/0"
  hap[is.na(hap)]<-"0/0"
  hap[grep("/",hap,invert=TRUE)]<-unlist(lapply(hap[grep("/",hap,invert = TRUE)],function(h){
    newh<-paste(h,h,sep="/")
  }))
  hap[grep("\\w+/\\w+/\\w+",hap)]<-"0/0" #remove any with multiple haplotypes (wtf?)
  gts<-as.data.frame(do.call(rbind,strsplit(hap,"/")))
  colnames(gts)<-c(idsa,idsb)
  return(gts)
}))
colnames(haps.keep)<-gsub("\\s+","",colnames(haps.keep))
rownames(haps.keep)<-colnames(haps.miss)[3:length(haps.miss)]
#remove consensus sequences
hkeep<-apply(haps.keep,2,function(x) !any(x=="consensus"))
haplotypes<-haps.keep[,hkeep==TRUE]
write.csv(haplotypes,"hap_genotypes.txt",row.names=TRUE,col.names=TRUE)
#large file
haplotypes<-read.csv("hap_genotypes.txt",row.names = 1,header = TRUE)

#remove any found in less than 90% of individuals
hap90<-prune.loci(haplotypes,0.1)
write.table(hap90,"HapsPolymorphicIn90.txt", quote=F,sep="\t",row.names=F)
#remove those loci not in hardy weinberg equilibrium
hwe<-data.frame()
for(x in seq(1,(ncol(hap90)-2),2)){
  hwe<-rbind(hwe,hwe.test(hap90[,c(x,(x+1))]))
}
keep.loci<-hwe[as.numeric(as.character(hwe$chi.result))>0.05,]
keep.names<-unlist(lapply(as.list(keep.loci$locus.name),function(x){
  names<-c(paste(x,"a",sep=""),paste(x,"b",sep=""))
  return(names) }))
hap.keep<-hap90[,colnames(hap90) %in% keep.names]

hap.keep$ID<-rownames(hap.keep)
hap.keep$sex<-gsub("sample_(\\w{3}).*","\\1",rownames(hap.keep))
hap.keep<-hap.keep[,c("ID","sex",keep.names)]
write.table(hap.keep,"dradPrunedHaps.txt",col.names=TRUE,row.names=TRUE,sep='\t',quote=FALSE)

#For CERVUS
s<-rep_cervus(nloci=2000,nreps=1,gens=hap.keep,out.prefix = "dradPrunedHaps") #for real analysis
c<-rep_cervus(nloci=c(50,100,150,300,200,400,800,1600),nreps=10,gens=hap.keep,out.prefix = "dradPrunedHaps")
write.table(hap.keep$ID[grep("FEM",hap.keep$ID)],
            "candidate_moms.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
off<-hap.keep$ID[grep("OFF",hap.keep$ID)]
dad<-hap.keep$ID[grep("PRM",hap.keep$ID)]
offs<-data.frame(offspring=off,KnownParent=NA,stringsAsFactors = FALSE)
for(i in 1:length(off)){
  d<-which(dadn==offn[i])
  if(length(d)==0){ offs[i,2]<-"" 
  }else { offs[i,2]<-as.character(dad[d]) }
}
write.table(offs,"offspring.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

#generate pairwise combos for band sharing
for(i in 1:length(rownames(hap.keep))){
  for(j in 1:length(rownames(hap.keep))){
    write.table(cbind(rownames(hap.keep)[c(i,j)]),"../relatedness/pairwise.combinations.txt",
                sep='\t',quote=FALSE,col.names = FALSE,row.names = FALSE,append = TRUE)
  }
}

## Create lists of males
#biallelic
bi<-read.delim("../parentage_biallelic/dradPruned200_1.txt")
write.table(grep("P",bi$ID,value = TRUE),"D:/SCA/parentage_biallelic/males.txt",
            quote=FALSE,col.names=FALSE,row.names = FALSE)
#haplotypes
hp<-read.delim("../parentage_haplotypes/dradPrunedHaps200_1.txt")
write.table(grep("P",hp$ID,value = TRUE),"D:/SCA/parentage/males.txt",
            quote=FALSE,col.names=FALSE,row.names = FALSE)

#######################PREVIOUS PRUNING
#prune to remove non-polymorphic loci
#<-apply(hapgen,2,function(x){ length(which(x=="consensus")) })
#all.poly<-consensus.num[consensus.num==0]
#poly<-hapgen[,colnames(hapgen) %in% names(all.poly)]

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