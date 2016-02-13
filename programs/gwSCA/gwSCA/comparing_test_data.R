#Author: Sarah P. Flanagan
#Date: 4 August 2015
#Purpose: Compare test input data to output from gwSCA

setwd("C:\\Users\\Sarah\\Documents\\Visual Studio 2013\\Projects\\gwSCA\\gwSCA")
in.alleles<-read.table("test_data.alleles.txt", header=T, sep='\t')
in.ind<-read.table("test_data.ind.txt", header=T)
out.alleles<-read.table("gw_sca_alleles.txt", header=T, sep='\t')
out.ind<-read.table("gw_sca_ind_info.txt", header=T)

 in.alleles[!(in.alleles[1,] %in% out.alleles[1,])]
 in.alleles[!(in.alleles[2,] %in% out.alleles[2,])]
 in.alleles[!(in.alleles[3,] %in% out.alleles[3,])]


