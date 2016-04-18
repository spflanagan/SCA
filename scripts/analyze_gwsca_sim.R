setwd("E:/ubuntushare/SCA/results/sca_simulation_output")
list.files(pattern="genotypes")
start<-read.table("ddraddist.test.afs.genotypes.start.txt")
end<-read.table("ddraddist.test.afs.genotypes.end.txt")

calc.maj.gf<-function(x){
	freq<-table(x)/sum(table(x))
	mgf<-max(freq)
	gen<-names(freq)[which.max(freq)]
	result<-rbind(gen,mgf)
	return(result)
}

start.gf<-apply(start[,-1],2,calc.maj.gf)
summary(t(start.gf))
# V1            V2      
# 0/0:8000   0.9724:8000  

end.gf<-apply(end[,-1],2,calc.maj.gf)
summary(t(end.gf))
h<-hist(as.numeric(end.gf[2,]),breaks=100)
plot(h$breaks[2:length(h$breaks)],h$density)
empirical.afs<-read.table("../biallelic/empirical_allelefreqs.txt")
lines(empirical.afs)