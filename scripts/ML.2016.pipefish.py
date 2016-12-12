## Genomic selection component analysis, John K. Kelly, June 2016
## Design based on field study of Iron Mountain population of M. guttatus conducted in 2014

## pipefish1:: models specific to data of Sarah Flannagan (July 2016)


def getDP(model,param_list):

	if model==0:
		q=param_list[0]
		x=[q*q,2*q*(1-q),1-q*q-2*q*(1-q)]#zygotes
		f0=param_list[1]
		f1=param_list[2]
		f2=param_list[3]
		wbar=f0*x[0]+f1*x[1]+f2*x[2]
		qm= f0*x[0]/wbar + 0.5*f1*x[1]/wbar
		DP=[qm-q]



def calc_v1(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped): 
	def scipy_ln_like3(x):
		return -ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,0, x,matplant,ungenotyped,1)

	bounds = [ (zy,1.0-zy), (minFec,maxFec), (LogSeedSD/2.0,2*LogSeedSD) ] # q,meanLnSeed,sig
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like3, x3ab, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like3(solution)
	solution.append(ln_l)
	zbob=ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant,ungenotyped,1)
	return solution

def calc_v0(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped): 
	def scipy_ln_like3(x):
		return -ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,0, x,matplant,ungenotyped,0)

	bounds = [ (zy,1.0-zy), (zy,1.0-zy), (minFec,maxFec), (LogSeedSD/2.0,2*LogSeedSD) ] # qM,qF,mu,sig
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like3, x3ab, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like3(solution)
	solution.append(ln_l)
	zbob=ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant,ungenotyped,0)
	return solution

def calc_v2(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped): 
	def scipy_ln_like3(x):
		return -ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,0, x,matplant,ungenotyped,2)

	bounds = [ (zy,1.0-zy), (minFec,maxFec),(minFec,maxFec),(minFec,maxFec), (LogSeedSD/2.0,2*LogSeedSD) ] #  q,mu1,mu2,mu3,sig
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like3, x3ab, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like3(solution)
	solution.append(ln_l)
	zbob=ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant,ungenotyped,2)
	return solution

def calc_v3(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped): 
	def scipy_ln_like3(x):
		return -ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,0, x,matplant,ungenotyped,3)

	bounds = [ (zy,1.0-zy), (zy,1.0-zy), (minFec,maxFec), (LogSeedSD/2.0,2*LogSeedSD) ] # q,qX,mu,sig
	best, val, d = optimize.fmin_l_bfgs_b(scipy_ln_like3, x3ab, approx_grad=True, bounds=bounds)
	solution = list(best)
	ln_l = -scipy_ln_like3(solution)
	solution.append(ln_l)
	zbob=ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,1, list(best),matplant,ungenotyped,3)
	return solution


def ln_like3(parents,famStr,RRL,RAL,AAL,FLnum,info, param_list,matplant,ungenotyped,model): 

# gfreqs[j] = pregnant male genotype frequencies
# FemFreqs[j] = wild caught female genotype frequencies
# mfreqs[j] = successful female genotype frequencies

	if model==1: # q,meanLnSeed,sig
		q=param_list[0]
		gRR=q*q
		gRA=2*q*(1-q)
		gAA=1-gRA-gRR
		mu0=param_list[1]
		mu1=param_list[1]
		mu2=param_list[1]
		sig=param_list[2]
		gfreqs=[gRR,gRA,gAA]
		mfreqs=[gRR,gRA,gAA]
		FemFreqs=[gRR,gRA,gAA]

	elif model==0: # qM,qF,mu,sig
		q=param_list[0] # qM
		gfreqs=[q*q,2*q*(1-q),(1-q)*(1-q)] 
		q=param_list[1] # qF
		FemFreqs=[q*q,2*q*(1-q),(1-q)*(1-q)]
		mfreqs=  [q*q,2*q*(1-q),(1-q)*(1-q)]
		mu0=param_list[2]
		mu1=param_list[2]
		mu2=param_list[2]
		sig=param_list[3]

	elif model==2: # q,mu1,mu2,mu3,sig
		q=param_list[0]
		gfreqs=  [q*q,2*q*(1-q),(1-q)*(1-q)] 
		FemFreqs=[q*q,2*q*(1-q),(1-q)*(1-q)]
		mfreqs=  [q*q,2*q*(1-q),(1-q)*(1-q)]
		mu0=param_list[1]
		mu1=param_list[2]
		mu2=param_list[3]
		sig=param_list[4]

	elif model==3: # q,qX,mu,sig
		q=param_list[0] # q
		gfreqs=[q*q,2*q*(1-q),(1-q)*(1-q)]
		FemFreqs=[q*q,2*q*(1-q),(1-q)*(1-q)]
		q=param_list[1] # qX
		mfreqs=  [q*q,2*q*(1-q),(1-q)*(1-q)]
		mu0=param_list[2]
		mu1=param_list[2]
		mu2=param_list[2]
		sig=param_list[3]


	ln_l = 0.0

	ll_seed=[0,0,0]
	for k in range(len(ungenotyped)): # lnl of seed counts for ungenotyped males
		ll_seed[0]= -log(sig**2.0)/2.0 - (ungenotyped[k]-mu0)**2.0 / (2*sig**2.0)
		ll_seed[1]= -log(sig**2.0)/2.0 - (ungenotyped[k]-mu1)**2.0 / (2*sig**2.0)
		ll_seed[2]= -log(sig**2.0)/2.0 - (ungenotyped[k]-mu2)**2.0 / (2*sig**2.0)
		Likelihood=0.0
		for j in range(3):
			Likelihood+= (gfreqs[j]*exp(ll_seed[j]))
		if Likelihood>0:
			ln_l+=log(Likelihood)
		else:
			return float('-inf')


	lnNM=ln_l # total for non-measured
	for j in range(parents):

		pdata=[RRL[j][0],RAL[j][0],AAL[j][0]] # prob all genotype data given MG 

		if FLnum[j][0]==0: # wild females
			Likelihood=0.0
			for k in range(3):
				Likelihood+= (FemFreqs[k]*pdata[k])
		elif FLnum[j][0]==1: # pregnant male
			scored_progeny=len(RRL[j])-1	
			if scored_progeny>0:
				sireSet=[] #how many distinct sires?
				for k in range(1,len(famStr[j])):
					if famStr[j][k] not in sireSet:
						sireSet.append(famStr[j][k])
				Product_over_sires=[1,1,1]	# index is mat genotype	
				#print j,len(sireSet),famStr[j]
				for zk in range(len(sireSet)):
				
					fmx=[[1.0,1.0,1.0],[1.0,1.0,1.0],[1.0,1.0,1.0]]
					for k in range(1,len(RRL[j])):
						if famStr[j][k] == sireSet[zk]: # in this set of full sibs

							fmx[0][0]*=(RRL[j][k])
							fmx[1][0]*=(0.5*RRL[j][k]+0.5*RAL[j][k]) 
							fmx[2][0]*= RAL[j][k]

							fmx[0][1]*=( 0.5*RRL[j][k]+0.5*RAL[j][k] )
							fmx[1][1]*=(0.25*RRL[j][k]+0.5*RAL[j][k]+0.25*AAL[j][k]) 
							fmx[2][1]*=( 0.5*RAL[j][k]+0.5*AAL[j][k] )

							fmx[0][2]*=( RAL[j][k] )
							fmx[1][2]*=(0.5*RAL[j][k]+0.5*AAL[j][k]) 
							fmx[2][2]*=( AAL[j][k] )  
														
					for matG in range(3):
						Product_over_sires[matG]*=(mfreqs[0]*fmx[matG][0]+mfreqs[1]*fmx[matG][1]+mfreqs[2]*fmx[matG][2])	

				for k in range(3):
					pdata[k]*=Product_over_sires[k]	


			ll_seed=[0,0,0]
			if FLnum[j][1]>=0.0:
				ll_seed[0]=-log(sig**2.0)/2.0 - (FLnum[j][1]-mu0)**2.0 / (2*sig**2.0)
				ll_seed[1]=-log(sig**2.0)/2.0 - (FLnum[j][1]-mu1)**2.0 / (2*sig**2.0)
				ll_seed[2]=-log(sig**2.0)/2.0 - (FLnum[j][1]-mu2)**2.0 / (2*sig**2.0)
				for k in range(3):
					if pdata[k]>0:
						pdata[k]= exp(log(pdata[k]) + ll_seed[k])
					else:
						pdata[k]=0.0

			Likelihood=0.0
			for k in range(3):
				Likelihood+= (gfreqs[k]*pdata[k])

		if Likelihood>0.0:
			ln_l+= log(Likelihood)
			if info==1 and VBS == 1:
				outb0.write(matplant[j]+'\t'+str(model)+'\t'+str(j)+'\t'+str(log(Likelihood))+'\n')
		else:
			return float('-inf')

	return ln_l



##### main program 
import sys
from scipy import optimize
from math import pow, log, exp
from scipy.stats.distributions import chi2

# GLOBAL CONSTANTS
zy=0.001
maximum_no_parents=500
#skip rare alleles
minPQ=0.05*(1.0-0.05)

VBS = 1 # set to 1 for detailed output
if VBS == 1:
	outb0 =open("Detail.output.txt","w")
	#outb1 =open("Detailed.1.txt","w")
	#outb2 =open("Detailed.2.txt","w")


gene = sys.argv[1] #sNNffold_1
maxf=0
in0  =open(sys.argv[1]+".LLm.txt", "rU")
out0 =open("ML.2016."+sys.argv[1]+".txt","w")
#out1 =open("S4.2014."+sys.argv[1]+".txt","w")


LnSeed={}
x=[0,0.0,0.0]
in1  =open("par.fit.pipefish.txt", "rU") # Log[seed+1] data
for line_idx, line in enumerate(in1):
	cols = line.replace('\n', '').split('\t')
	# 338	2.8082109729
	LnSeed[cols[0]]=float(cols[1])
	x[0]+=1
	x[1]+=float(cols[1])
	x[2]+=float(cols[1])**2.0
mean=x[1]/x[0]
maxFec=5.0*mean
minFec=0.01*mean
var= x[2]/x[0] - mean**2.0
var= var*float(x[0])/float(x[0]-1)
LogSeedSD=var**(0.5)
#print "Mean and SD of seed ",mean,LogSeedSD

for line_idx, line in enumerate(in0):
	cols = line.replace('\n', '').split('\t')
	# LG1	4222203	FEM.001	A,A	1	0.0001	0.0001

	if line_idx==0:
		matplant=[]
		famStr =[[] for j in range(maximum_no_parents)]
		RRL =[[] for j in range(maximum_no_parents)]
		RAL =[[] for j in range(maximum_no_parents)]
		AAL =[[] for j in range(maximum_no_parents)]
		FLnum=[[-9,-9] for j in range(maximum_no_parents)]
		parents=0
		total_plants=0
		found_not=[0,0]
		snpX=cols[1]

	if cols[1] != snpX: # new snp		
		# make list of fecundities for ungenotyped pregnant males
		ungenotyped=[]
		for zz in LnSeed:
			if zz not in matplant:
				ungenotyped.append(LnSeed[zz])
		if VBS == 1:
			outb0.write(snpX+'\n')
		# Analzye data
		x3ab=[0.5,mean,LogSeedSD] # q,meanLnSeed,sig
		global_v1 = calc_v1(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) # no selection
		px=global_v1[0]
		print "Pregnant males ",found_not," ll0 ",global_v1[3]
		if px*(1-px)<minPQ:
			pass # skip rare alleles
		else:
			out0.write(cols[0]+'\t'+snpX+'\t'+str(parents))
			for yy in range(len(global_v1)):
				out0.write( '\t'+str(round(global_v1[yy],4)) )
 
			x3ab=[px,px,global_v1[1],global_v1[2]] # qM,qF,mu,sig 
			global_v0 = calc_v0(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
			global_v0.append(2.0*(global_v0[4]-global_v1[3]))
			for yy in range(len(global_v0)):
				out0.write( '\t'+str(round(global_v0[yy],4)) )

			x3ab=[px,global_v1[1],global_v1[1],global_v1[1],global_v1[2]] # q,mu1,mu2,mu3,sig 
			global_v2 = calc_v2(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
			global_v2.append( 2.0*(global_v2[5]-global_v1[3]) )
			for yy in range(len(global_v2)):
				out0.write( '\t'+str(round(global_v2[yy],4)) )

			x3ab=[px,px,global_v1[1],global_v1[2]] # q,qX,mu,sig 
			global_v3 = calc_v3(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
			global_v3.append(2.0*(global_v3[4]-global_v1[3]))
			for yy in range(len(global_v3)):
				out0.write( '\t'+str(round(global_v3[yy],4)) )
				
			out0.write('\n')

		#reset variables for data entry
		matplant=[]
		famStr =[[] for j in range(maximum_no_parents)]
		RRL =[[] for j in range(maximum_no_parents)]
		RAL =[[] for j in range(maximum_no_parents)]
		AAL =[[] for j in range(maximum_no_parents)]
		FLnum=[[-9,-9] for j in range(maximum_no_parents)]
		parents=0
		total_plants=0
		found_not=[0,0]
		snpX=cols[1]
	# sNNffold_1	17764	loser.1-4-1	G,T	1.0	0.5	1e-06
	# sNNffold_1	17764	1-3-2	G,T	2	1.0	1.0	1.0	2	1.0	0.34	1e-06	0	1.0	0.34	1e-06

	matplant.append(cols[2])
	if cols[2].split(".")[0]=="FEM":
		FLnum[parents][0]=0

		LG=[float(cols[4]),float(cols[5]),float(cols[6])]
		if LG[0]+LG[1]+LG[2]>0:
			RRL[parents].append(LG[0])
			RAL[parents].append(LG[1])
			AAL[parents].append(LG[2])
		else:
			RRL[parents].append(1.0)
			RAL[parents].append(1.0)
			AAL[parents].append(1.0)

	else:
		FLnum[parents][0]=1
		try:
			FLnum[parents][1]=LnSeed[cols[2]] # seed set 
			found_not[0]+=1
		except KeyError:
			found_not[1]+=1
			#pass

		scoredProgeny=int(cols[4])
		LG=[float(cols[5]),float(cols[6]),float(cols[7])]	# parent likelihoods
		famStr[parents].append(-9) # place holder
		if LG[0]+LG[1]+LG[2]>0: # not sure this is needed
			RRL[parents].append(LG[0])
			RAL[parents].append(LG[1])
			AAL[parents].append(LG[2])
		else:
			RRL[parents].append(1.0)
			RAL[parents].append(1.0)
			AAL[parents].append(1.0)

		for j in range(scoredProgeny):
			famStr[parents].append(cols[8+4*j]) # sire IDs
			LG=[float(cols[9+4*j]),float(cols[10+4*j]),float(cols[11+4*j])]
			if LG[0]+LG[1]+LG[2]>0:
				RRL[parents].append(LG[0])
				RAL[parents].append(LG[1])
				AAL[parents].append(LG[2])
			else:
				RRL[parents].append(1.0)
				RAL[parents].append(1.0)
				AAL[parents].append(1.0)

	parents+=1


# last snp
ungenotyped=[]
for zz in LnSeed:
	if zz not in matplant:
		ungenotyped.append(LnSeed[zz])

if VBS == 1:
	outb0.write(snpX+'\n')
# Analzye data
x3ab=[0.5,mean,LogSeedSD] # q,meanLnSeed,sig
global_v1 = calc_v1(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) # no selection
px=global_v1[0]
if px*(1-px)<minPQ:
	pass # skip rare alleles
else:
	out0.write(cols[0]+'\t'+snpX+'\t'+str(parents))
	for yy in range(len(global_v1)):
		out0.write( '\t'+str(round(global_v1[yy],4)) )

	x3ab=[px,px,global_v1[1],global_v1[2]] # qM,qF,mu,sig 
	global_v0 = calc_v0(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
	global_v0.append(2.0*(global_v0[4]-global_v1[3]))
	for yy in range(len(global_v0)):
		out0.write( '\t'+str(round(global_v0[yy],4)) )

	x3ab=[px,global_v1[1],global_v1[1],global_v1[1],global_v1[2]] # q,mu1,mu2,mu3,sig 
	global_v2 = calc_v2(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
	global_v2.append( 2.0*(global_v2[5]-global_v1[3]) )
	for yy in range(len(global_v2)):
		out0.write( '\t'+str(round(global_v2[yy],4)) )

	x3ab=[px,px,global_v1[1],global_v1[2]] # q,qX,mu,sig 
	global_v3 = calc_v3(parents,x3ab,famStr,RRL,RAL,AAL,FLnum,matplant,ungenotyped) 
	global_v3.append(2.0*(global_v3[4]-global_v1[3]))
	for yy in range(len(global_v3)):
		out0.write( '\t'+str(round(global_v3[yy],4)) )
		
	out0.write('\n')

			
