#pragma once

#include "random_numbers.h"
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <vector>
#include <array>
#include <sstream>

using namespace std;

class ld_info
{
public:
	double d;
	double dprime;
};

class chromosome_emulator
{
public:
	vector<double> locus_emulator;

	chromosome_emulator()
	{
		locus_emulator = vector<double>();
	}
};

class allele_freqs
{
public:
	vector<vector<double>> freqs;

	allele_freqs()
	{
		freqs = vector<vector<double>>();
	}
};

class locus_statistics
{
public:
	vector<chromosome_emulator> fst;
	vector<chromosome_emulator> exp_het;
	vector<chromosome_emulator> weighted_fst;
	vector<chromosome_emulator> locus_id;
	vector<allele_freqs> pop1_freqs;
	vector<allele_freqs> pop2_freqs;
	vector<vector<int>> num_alleles;

	locus_statistics()
	{
		fst = vector<chromosome_emulator>();
		exp_het = vector<chromosome_emulator>();
		weighted_fst = vector<chromosome_emulator>();
		locus_id = vector<chromosome_emulator>();
		num_alleles = vector<vector<int>>();
		pop1_freqs = vector<allele_freqs>();
		pop2_freqs = vector<allele_freqs>();
	}

	void stats_initialize(int num_chrom, int num_markers, int alleles)
	{
		int j, jj, jjj;
		for (j = 0; j < num_chrom; j++)
		{ 
			exp_het.push_back(chromosome_emulator());
			fst.push_back(chromosome_emulator());
			locus_id.push_back(chromosome_emulator());
			num_alleles.push_back(vector<int>());
			pop1_freqs.push_back(allele_freqs());
			pop2_freqs.push_back(allele_freqs());
			for (jj = 0; jj < num_markers; jj++)
			{
				exp_het[j].locus_emulator.push_back(double());
				fst[j].locus_emulator.push_back(double());
				locus_id[j].locus_emulator.push_back(double());
				num_alleles[j].push_back(int());
				pop1_freqs[j].freqs.push_back(vector<double>());
				pop2_freqs[j].freqs.push_back(vector<double>());
				for (jjj = 0; jjj < alleles; jjj++)
				{
					pop1_freqs[j].freqs[jj].push_back(0);
					pop2_freqs[j].freqs[jj].push_back(0);
				}
			}
		}
	}
};

class chromosome
{
public:
	vector<int> loci;
	vector<double> allelic_effects;

	chromosome()
	{
		loci = vector<int>();
		allelic_effects = vector<double>();
	}
};//chromosome

class individual
{
public:
	vector<chromosome> maternal;
	vector<chromosome> paternal;
	double genotype, phenotype;
	bool alive, female;
	int mate;
	vector<int> offspring_index;

	individual()
	{
		maternal = vector<chromosome>();
		paternal = vector<chromosome>();
		genotype = double();
		phenotype = double();
		alive = bool();
		female = bool();
		mate = int();
		offspring_index = vector<int>();
	}

	void initialize_individual(int numchrom, int nummarkers, int numqtl)
	{
		int j, jj;
		alive = true;
		if (genrand() < 0.5)
			female = true;
		else
			female = false;
		for (j = 0; j < numchrom; j++)
		{
			maternal.push_back(chromosome());
			paternal.push_back(chromosome());
			for (jj = 0; jj < nummarkers; jj++)
			{
				maternal[j].loci.push_back(int());
				paternal[j].loci.push_back(int());
			}
			for (jj = 0; jj < numqtl; jj++)
			{
				maternal[j].allelic_effects.push_back(double());
				paternal[j].allelic_effects.push_back(double());
			}
		}
	}
};//end individual

class population
{
public:
	//declare all variables
	vector <individual> adults;
	vector<individual> offspring;
	int pop_size, max_fecundity, max_encounters, max_num_prog;
	int num_mal, num_fem, num_off;
	bool extinct, bootstrap, prior_qtl;
	double gaussian_preference_mean;
	double significance_level, sigma;
	double mean_mal_trait, mean_fem_trait;
	double recombination_rate, avg_pairwise_d, avg_ld;
	int peak_window, num_ld;
	int mutation_rate, mutational_variance;
	int environmental_variance, environmental_sd, allelic_sd;
	ifstream afs_file;
	string afs_file_name;

	//variable parameters
	int carrying_cap, num_gen;
	int num_chrom, num_markers, num_qtl, num_alleles;
	int adult_samplesize, offspring_samplesize, male_samplesize, female_samplesize, dad_samplesize;

	//things that help keep track of other things
	vector<chromosome_emulator> qtl_tracker;
	vector<int> qtl_list;
	vector<chromosome_emulator> primary_frequency, major_allele;
	vector<int> sampled_adults, sampled_males, sampled_females, sampled_off, sampled_dads;
	vector<ld_info> pop_ld;
	vector<vector<chromosome>> inferred_dads;
	locus_statistics adult_offspring, male_female, males_dads;


	void set_parameters()
	{
		carrying_cap = 5000;
		num_gen = 200;

		num_chrom = 4;
		num_markers = 2000;
		num_qtl = 2;
		num_alleles = 2;
		prior_qtl = false;
		afs_file_name = "E://ubuntushare//SCA//results//biallelic//empirical_allelefreqs.txt";

		adult_samplesize = 440;
		offspring_samplesize = 160;
		male_samplesize = 87;//females in pipefish
		female_samplesize = 196;//males in pipefish
		dad_samplesize = 160;//moms in pipefish

		recombination_rate = 0.2;
		environmental_variance = 0;
		environmental_sd = sqrt(environmental_variance);
		mutation_rate = 0;
		mutational_variance = 0;
		max_fecundity = 4;
		max_num_prog = max_fecundity*carrying_cap;
		max_encounters = 50;
		gaussian_preference_mean = 4;
		allelic_sd = 0.5;
		
		significance_level = 0.5;
		peak_window = 50;
		num_ld = 2000;
		sigma = 15;
		bootstrap = false;

		extinct = false;
		sgenrand(time(0));
	}//end set_parameters

	void use_empirical_afs(vector<double>& tempallele)
	{
		string line;
		vector<double> intervals;
		vector<double> densities; 
		vector<int> count_intervals;
		double t, tt, last_interval, last_prob;
		int j, jj, jjj, int_index, count;

		afs_file.open(afs_file_name);
		while (getline(afs_file, line))
		{
			stringstream temp;
			temp.str(line);
			temp >> t >> tt;
			intervals.push_back(t);
			densities.push_back(tt);
		}
		
		for (j = 0; j < num_alleles; j++)
		{
			t = genrand();
			int_index = 0;
			for (jj = 0; jj < intervals.size(); jj++)
			{				
				if (intervals[jj] == t)
					int_index = jj;	
			}
			tt = randnorm(densities[int_index], intervals[0] / 2);
			tempallele.push_back(tt);
		}

		vector<int> alleles;
		vector<double> allele_probs;
		for (j = 0; j < num_alleles; j++)
			alleles.push_back(carrying_cap%num_alleles);
			
		int num_allele_intervals = intervals.size() / num_alleles;
		if (num_allele_intervals*num_alleles < intervals.size())
			num_allele_intervals++;
		for (j = 0; j < num_allele_intervals; j++)
			allele_probs.push_back(0);
		int_index = 0;
		count_intervals.push_back(0);
		for (j = 0; j < intervals.size(); j++)
		{
			if (int_index >= num_allele_intervals)
			{
				count_intervals.push_back(0);
				int_index++;
			}
			allele_probs[int_index] = allele_probs[int_index]+ densities[j];
			count_intervals[int_index]++;
		}
		for (j = 0; j < num_allele_intervals; j++)
			allele_probs[j] = allele_probs[j] / count_intervals[j];
		for (j = 0; j < carrying_cap; j++)
		{
			for (jj = 0; jj < num_chrom; jj++)
			{
				for (jjj = 0; jjj < num_markers; jjj++)
				{
					t = genrand();
					int_index = count = 0;
					for (jj = 0; jj < intervals.size(); jj++)
					{
						if (count >= num_allele_intervals)
							count++;
						if (intervals[jj] == t)
							int_index = count;
					}
					tt = randnorm(allele_probs[int_index], intervals[0] / 2);
					adults[j].maternal[jj].loci[jjj] = j%num_alleles;
					adults[j].paternal[jj].loci[jjj] = j%num_alleles;
				}
			}
		}

	}
	
	

	void initialize(bool known_qtls, bool empirical_afs)
	{
		cout << "Initializing the population.\n";
		int j, jj, jjj, count;
		num_fem = num_mal = 0;
		mean_fem_trait = mean_mal_trait = 0;
		for (j = 0; j < carrying_cap; j++)
		{
			adults.push_back(individual());
			adults[j].initialize_individual(num_chrom, num_markers, num_alleles);
		}
		
		//tracker structures
		adult_offspring.stats_initialize(num_chrom, num_markers, num_alleles);
		male_female.stats_initialize(num_chrom, num_markers, num_alleles);
		males_dads.stats_initialize(num_chrom, num_markers, num_alleles);
		for (j = 0; j < num_chrom; j++)
		{
			qtl_tracker.push_back(chromosome_emulator());
			primary_frequency.push_back(chromosome_emulator());
			major_allele.push_back(chromosome_emulator());
			pop_ld.push_back(ld_info());
			for (jj = 0; jj < num_markers; jj++)
			{
				primary_frequency[j].locus_emulator.push_back(double());
				major_allele[j].locus_emulator.push_back(double());
				qtl_tracker[j].locus_emulator.push_back(double());
			}
		}
		for (j = 0; j < dad_samplesize; j++)
		{
			inferred_dads.push_back(vector<chromosome>());
			for (jj = 0; jj < num_chrom; jj++)
			{
				inferred_dads[j].push_back(chromosome());
				for (jjj = 0; jjj < num_markers; jjj++)
					inferred_dads[j][jj].loci.push_back(int());
				for (jjj = 0; jjj < num_qtl; jjj++)
					inferred_dads[j][jj].allelic_effects.push_back(double());
			}
		}

		//assign genotypes
		if (known_qtls)
		{
			int known_qtls[] = { 444, 82, 393, 218, 699, 477, 603, 694 };
			int total_num_qtls = (num_qtl*num_chrom);
			qtl_list.assign(known_qtls, known_qtls+total_num_qtls);
		}
		else
		{
			for (j = 0; j < (num_qtl*num_chrom); j++)
				qtl_list.push_back(randnum(num_markers));
		}
		count = 0;
		for (j = 0; j < num_chrom; j++)
		{
			for (jj = 0; jj < num_markers; jj++)
			{
				if (qtl_list[count] == jj)
					qtl_tracker[j].locus_emulator[jj] = 1;
				else
					qtl_tracker[j].locus_emulator[jj] = 0;
			}
		}//end qtls
		
		vector<double> tempallele;
		vector<double> tempmarkers;
		int index;
		for (j = 0; j < num_alleles; j++)
			tempallele.push_back(randnorm(0, allelic_sd));
		if (empirical_afs)
		{
			for (j = 0; j < num_alleles; j++)
				tempmarkers.push_back(j%num_alleles);
		}
		for (j = 0; j < carrying_cap; j++)
		{
			adults[j].phenotype = adults[j].genotype = 0;
			if (empirical_afs)
				index = d_empirical_afs(num_alleles, afs_file_name);	
			for (jj = 0; jj < num_chrom; jj++)
			{
				
				for (jjj = 0; jjj < num_qtl; jjj++)
				{
					if (!empirical_afs)
					{
						adults[j].maternal[jj].allelic_effects[jjj] = tempallele[j%num_alleles];
						adults[j].paternal[jj].allelic_effects[jjj] = tempallele[j%num_alleles];
					}
					else
					{
						adults[j].maternal[jj].allelic_effects[jjj] = tempallele[index];
						adults[j].paternal[jj].allelic_effects[jjj] = tempallele[index];
					}
					adults[j].genotype = adults[j].genotype + adults[j].maternal[jj].allelic_effects[jjj]
						+ adults[j].paternal[jj].allelic_effects[jjj];
				}
				adults[j].phenotype = adults[j].genotype + randnorm(0, environmental_sd);
				for (jjj = 0; jjj < num_markers; jjj++)
				{
					if (!empirical_afs)
					{
						adults[j].maternal[jj].loci[jjj] = j%num_alleles;
						adults[j].paternal[jj].loci[jjj] = j%num_alleles;
					}
					else
					{
						adults[j].maternal[jj].loci[jjj] = tempmarkers[index];
						adults[j].paternal[jj].loci[jjj] = tempmarkers[index];
					}
				}

			}
			if (adults[j].female == true)
			{
				mean_fem_trait = mean_fem_trait + adults[j].phenotype;
				num_fem++;
			}
			else
			{
				mean_mal_trait = mean_mal_trait + adults[j].phenotype;
				num_mal++;
			}
		}
		pop_size = carrying_cap;
		mean_fem_trait = mean_fem_trait / num_fem;
		mean_mal_trait = mean_mal_trait / num_mal;

		
	}//end initialize

	int determine_pop_size()
	{
		int count, p;
		count = 0;
		for (p = 0; p < carrying_cap; p++)
		{
			if (adults[p].alive)
				count++;
		}
		return count;
	}//end DeterminePopSize

	void RecombineChromosome(chromosome &RecombinedChr, individual &Parent, int WhichChromosome, double ExpectedRecombEvents)
	{
		int RCi, RCj;
		int NumberRecombEvents = 0;
		vector<int> SegmentStart, SegmentEnd;
		vector<int> BreakPoint;

		if (ExpectedRecombEvents < 6)
			NumberRecombEvents = poissonrand(ExpectedRecombEvents);
		if (ExpectedRecombEvents > 6)
			NumberRecombEvents = positiveroundnorm(ExpectedRecombEvents, sqrt(ExpectedRecombEvents));
		vector<bool> SegmentMaternal;
		int NumberSegments;
		bool StartMaternal;
		if (NumberRecombEvents > 20)
			NumberRecombEvents = 20;

		if (NumberRecombEvents > 0)
		{
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
				BreakPoint.push_back(randnum(num_markers));
			//sort breakpoints
			sort(begin(BreakPoint), end(BreakPoint));

			//first segment maternal or paternal?
			if (genrand() < 0.5)
				StartMaternal = true;
			else
				StartMaternal = false;

			NumberSegments = 1;
			SegmentStart.push_back(0);
			SegmentMaternal.push_back(StartMaternal);
			for (RCi = 0; RCi < NumberRecombEvents; RCi++)
			{
				SegmentEnd.push_back(BreakPoint[RCi]);
				SegmentStart.push_back(BreakPoint[RCi]);
				if (SegmentMaternal[RCi])
					SegmentMaternal.push_back(false);
				else
					SegmentMaternal.push_back(true);
				NumberSegments++;
			}//end RCi
			SegmentEnd.push_back(num_markers);

			//now pass allelic info to recombined chromosome
			for (RCi = 0; RCi < NumberSegments; RCi++)
			{
				if (SegmentMaternal[RCi])
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
						RecombinedChr.loci[RCj] = Parent.maternal[WhichChromosome].loci[RCj];
				}
				else
				{
					for (RCj = SegmentStart[RCi]; RCj < SegmentEnd[RCi]; RCj++)
						RecombinedChr.loci[RCj] = Parent.paternal[WhichChromosome].loci[RCj];
				}
			}//end RCi

			//now do the QTLs
			int count = 0;
			for (RCi = 0; RCi < num_qtl; RCi++)
			{
				for (RCj = 0; RCj < NumberSegments; RCj++)
				{
					if (qtl_list[count] >= SegmentStart[RCj] && qtl_list[count] < SegmentEnd[RCj])
					{
						if (SegmentMaternal[RCj])
							RecombinedChr.allelic_effects[RCi] = Parent.maternal[WhichChromosome].allelic_effects[RCi];
						else
							RecombinedChr.allelic_effects[RCi] = Parent.paternal[WhichChromosome].allelic_effects[RCi];
					}
				}
				count++;
			}//end RCi
		}//end numb recomb events > 0
		else
		{
			//No recombination
			if (genrand() < 0.5)
			{
				for (RCi = 0; RCi < num_markers; RCi++)
					RecombinedChr.loci[RCi] = Parent.maternal[WhichChromosome].loci[RCi];
				for (RCi = 0; RCi < num_qtl; RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.maternal[WhichChromosome].allelic_effects[RCi];
			}
			else
			{
				for (RCi = 0; RCi < num_markers; RCi++)
					RecombinedChr.loci[RCi] = Parent.paternal[WhichChromosome].loci[RCi];
				for (RCi = 0; RCi < num_qtl; RCi++)
					RecombinedChr.allelic_effects[RCi] = Parent.paternal[WhichChromosome].allelic_effects[RCi];
			}
		}//else (no recomb)
	}//end RecombineChromosome
	
	void mating(double GaussianPrefVariance)
	{
		int cc;
		int NumProg, prog_fem, prog_mal;
		double dubrand;
		NumProg = 0;
		//Mate Choice:
		int mm, nn, males;
		int Encounters, MaleID, irndnum, Females;
		double MeanFemaleTrait;
		bool MateFound, MalesPresent;
		double MeanMaleTrait, MateProb;
		int MaleIndex, FemaleID;
		int counter1, counter2, counter3;
		vector <int> MaleList;
		Females = 0;
		prog_fem = prog_mal = 0;
		num_mal = num_fem = 0;
		MalesPresent = false;
		MeanMaleTrait = 0;
		double SDMaleTrait = 0;
		MeanFemaleTrait = 0;
		counter1 = counter2 = counter3 = 0;
		//determine mean male trait
		offspring.resize(0);
		for (males = 0; males < pop_size; males++)
		{
			adults[males].mate = 0;
			adults[males].offspring_index.resize(0);
			if (!adults[males].female)
			{
				MalesPresent = true;
				MaleList.push_back(males);
				num_mal++;
				MeanMaleTrait = MeanMaleTrait + adults[males].phenotype;
			}
		} // end of males
		if (num_mal > 0) 
		{
			MeanMaleTrait = MeanMaleTrait / num_mal;
		}
		else 
		{
			MeanMaleTrait = 0;
			extinct = true;
		}
		for (males = 0; males < pop_size; males++)
		{
			if (!adults[males].female)
				SDMaleTrait = SDMaleTrait + (adults[males].phenotype - MeanMaleTrait)*(adults[males].phenotype - MeanMaleTrait);
			else
				adults[males].offspring_index.resize(0);
		}
		SDMaleTrait = sqrt(SDMaleTrait / num_mal);
		for (mm = 0; mm < pop_size; mm++)
		{
			MateFound = false;
			if (adults[mm].female && MalesPresent)
			{
				FemaleID = mm;
				MeanFemaleTrait = MeanFemaleTrait + adults[mm].phenotype;
				Females++;
				Encounters = 0;
				while (!MateFound && Encounters <= max_encounters)
				{
					irndnum = randnum(num_mal);
					MaleIndex = MaleList[irndnum];
					if (GaussianPrefVariance > 0)
						MateProb = exp(-0.5 * (adults[MaleIndex].phenotype - gaussian_preference_mean)*
							(adults[MaleIndex].phenotype - gaussian_preference_mean) / GaussianPrefVariance);
					else
						MateProb = 1;
					dubrand = genrand();
					if (dubrand < MateProb)
					{
						MateFound = true;
						MaleID = MaleIndex;
						adults[MaleID].mate++;
						counter2++;
						//Mated++;
					}
					Encounters++;
				}//while
				if (MateFound)
				{
					counter1++;
					adults[mm].mate++;
					if (NumProg >= max_num_prog)
						NumProg = max_num_prog - 1;
					//mother is Parent 1, mm
					//father is parent 2, MateID
					for (nn = 0; nn < max_fecundity; nn++)
					{
						adults[FemaleID].offspring_index.push_back(NumProg);
						adults[MaleID].offspring_index.push_back(NumProg);
						offspring.push_back(individual());
						offspring[NumProg].initialize_individual(num_chrom, num_markers, num_qtl);
						//calculate phenotype in mutation once the genotype is determined
						for (cc = 0; cc < num_chrom; cc++)//go through chromosome by chromosome
						{
							RecombineChromosome(offspring[NumProg].maternal[cc], adults[FemaleID], cc, recombination_rate);
							RecombineChromosome(offspring[NumProg].paternal[cc], adults[MaleID], cc, recombination_rate);						
						}//end of chromosome
						if (genrand() < 0.5)
						{
							offspring[NumProg].female = true;
							prog_fem++;
						}
						else
						{
							offspring[NumProg].female = false;
							prog_mal++;
						}
						NumProg++;
					}//for nn
				}
				if (!MateFound)
				{//Keep track of the females that didn't mate
					counter3++;
				}
			}//if female
		}//end of mm
		num_off = NumProg;
		num_fem = Females;
		MeanFemaleTrait = MeanFemaleTrait / Females;
	}//end mating

	void mutation()
	{
		int m, mm, irand, irand2, gg, ggg, irand3, locus;
		double rnd1, rnd2;
		double IndMutationRate;
		double MutSD;
		bool mutated;

		MutSD = sqrt(mutational_variance);
		IndMutationRate = mutation_rate * 2 * num_markers * num_chrom;

		for (m = 0; m < num_off; m++)
		{
			rnd1 = genrand();
			offspring[m].phenotype = 0;
			offspring[m].genotype = 0;
			mutated = false;
			if (rnd1 < IndMutationRate)
			{
				irand = randnum(num_chrom);
				irand2 = randnum(num_markers);
				locus = irand2;
				if (genrand() < 0.5)//affects maternal chromosome	
				{
					while (!mutated)
					{
						irand3 = randnum(num_alleles);
						if (!offspring[m].maternal[irand].loci[irand2] == irand3)
						{
							offspring[m].maternal[irand].loci[irand2] = irand3;
							mutated = true;
						}
					}
					if (qtl_tracker[irand].locus_emulator[irand2] == 1)
					{
						for (mm = irand*num_qtl; mm < num_qtl; mm++)
							offspring[m].maternal[irand].allelic_effects[mm] =
							offspring[m].maternal[irand].allelic_effects[mm] + randnorm(0, MutSD);
					}
				}
				else//affects paternal chromosome
				{
					while (!mutated)
					{
						irand3 = randnum(num_alleles);
						if (!offspring[m].paternal[irand].loci[locus] == irand3)
						{
							offspring[m].paternal[irand].loci[locus] = irand3;
							mutated = true;
						}
					}
					if (qtl_tracker[irand].locus_emulator[irand2] == 1)
					{
						for (mm = irand*num_qtl; mm < num_qtl; mm++)
							offspring[m].paternal[irand].allelic_effects[mm] =
							offspring[m].paternal[irand].allelic_effects[mm] + randnorm(0, MutSD);
					}
				}//end of if

				for (gg = 0; gg < num_chrom; gg++)
				{
					for (ggg = 0; ggg < num_qtl; ggg++){
						offspring[m].genotype = offspring[m].genotype +
							offspring[m].maternal[gg].allelic_effects[ggg] + offspring[m].paternal[gg].allelic_effects[ggg];
					}
				}
				offspring[m].phenotype = offspring[m].genotype + randnorm(0, environmental_sd);
			}
		}//end of m
	}//mutation

	void viability_selection(double selection_strength)
	{
		int j, prog_alive, mal_count;
		double survive_probability, drnum1, optimum, phen_sd, phen_mean, num;
		mal_count = prog_alive = 0;
		phen_mean = phen_sd = 0;
		//calc mean male phenotype and std dev
		for (j = 0; j < num_off; j++)
		{
			if (offspring[j].female)
			{
				phen_mean = phen_mean + offspring[j].phenotype;
				mal_count++;
			}
		}
		num = mal_count;
		phen_mean = phen_mean / num;
		for (j = 0; j < num_off; j++)
		{
			if (!offspring[j].female)
				phen_sd = phen_sd + (offspring[j].phenotype - phen_mean)*(offspring[j].phenotype - phen_mean);
		}
		phen_sd = sqrt(phen_sd / num);

		optimum = 0;
		for (j = 0; j < num_off; j++)
		{
			if (offspring[j].female)
			{
				offspring[j].alive = true;
				prog_alive++;
			}
			else
			{
				if (selection_strength > 0)
					survive_probability = exp(-1 * (offspring[j].phenotype - optimum)*(offspring[j].phenotype - optimum)
					/ (2 * selection_strength));
				else
					survive_probability = 1;
				drnum1 = genrand();
				if (drnum1 < survive_probability)
				{
					offspring[j].alive = true;
					prog_alive++;
				}
				else
					offspring[j].alive = false;
			}
		}
	}//end viability selection

	void density_regulation()
	{
		int p, pp, ppp;
		int iNumAdultsChosen;
		double CarCapUnfilled, ProgLeft, KeepProb;
		double DRrandnum;
		num_mal = 0;
		num_fem = 0;
		ProgLeft = 0;
		//count the ones that are still alive
		for (p = 0; p < num_off; p++)
		{
			if (offspring[p].alive)
				ProgLeft++;
		}
		CarCapUnfilled = carrying_cap;
		iNumAdultsChosen = 0;
		for (p = 0; p < num_off; p++)
		{
			if (offspring[p].alive)
			{
				if (ProgLeft == 0)
					KeepProb = 0;
				else
					KeepProb = CarCapUnfilled / ProgLeft;
				DRrandnum = genrand();
				if (DRrandnum<KeepProb)
				{//then turn it into an adult
					adults[iNumAdultsChosen].alive = true;
					adults[iNumAdultsChosen].mate = 0;
					for (pp = 0; pp < num_chrom; pp++)
					{
						for (ppp = 0; ppp < num_markers; ppp++)
						{
							adults[iNumAdultsChosen].maternal[pp].loci[ppp] = offspring[p].maternal[pp].loci[ppp];
							adults[iNumAdultsChosen].paternal[pp].loci[ppp] = offspring[p].paternal[pp].loci[ppp];
						}
						for (ppp = 0; ppp < num_qtl; ppp++)
						{
							adults[iNumAdultsChosen].maternal[pp].allelic_effects[ppp] = offspring[p].maternal[pp].allelic_effects[ppp];
							adults[iNumAdultsChosen].paternal[pp].allelic_effects[ppp] = offspring[p].paternal[pp].allelic_effects[ppp];
						}
					}
					adults[iNumAdultsChosen].phenotype = offspring[p].phenotype;
					adults[iNumAdultsChosen].genotype = offspring[p].genotype;
					if (offspring[p].female){
						adults[iNumAdultsChosen].female = true;
						num_fem++;
					}
					else{
						adults[iNumAdultsChosen].female = false;
						num_mal++;
					}
					CarCapUnfilled = CarCapUnfilled - 1;
					iNumAdultsChosen++;
				}//end of if KeepProb
				else
					offspring[p].alive = false;
			}//end of if Alive
			ProgLeft = ProgLeft - 1;
		}//end of for p
		pop_size = iNumAdultsChosen;
		if (pop_size== 0)
			extinct = true;
	}//end Density Regulation

	void calc_mean_trait_values()
	{
		int p;
		num_mal = num_fem = mean_mal_trait = mean_fem_trait = 0;
		for (p = 0; p < pop_size; p++)
		{
			if (adults[p].female)
			{
				mean_fem_trait = mean_fem_trait+ adults[p].phenotype;
				num_fem++;
			}
			if (!adults[p].female)
			{
				mean_mal_trait = mean_mal_trait + adults[p].phenotype;
				num_mal++;
			}
		}
		mean_mal_trait = mean_mal_trait / double(num_mal);
		mean_fem_trait = mean_fem_trait/ double(num_fem);
	}//end calc mean trait values

	ld_info adult_pop_ld(int chromosome_1, int chromosome_2, int allele_1, int allele_2)
	{
		//To do this, need to compare allele frequencies. 
		//D=x11-p1q1
		//x11 = observed freq of major allele in locus A and major allele in locus B
		//p1 = observed freq of major allele in locus A
		//q1 = observed freq of major allele in locus B
		//if D > 0, Dmax = min(p1q1, p2q2)
		//if D < 0, Dmax = min(p1q2, p2q1)
		//D' = D/Dmax
		ld_info result;
		double dNadult = pop_size;
		vector <int> num_allele_1;
		vector <int> num_allele_2;
		vector <double> freq_allele_1;
		vector <double> freq_allele_2;
		double **joint_alleles = new double *[num_alleles];
		double **D = new double *[num_alleles];
		double **Dmax = new double *[num_alleles];
		int count_a = 0;
		int f, ff, fff, al, count;
		double Dprime, d_allele_avgs;
		int numA1B1 = 0;

		for (al = 0; al < num_alleles; al++)
		{
			joint_alleles[al] = new double[num_alleles];
			D[al] = new double[num_alleles];
			Dmax[al] = new double[num_alleles];
		}
		//figure out which loci are polymorphic
		for (al = 0; al < num_alleles; al++)
		{
			num_allele_1.push_back(0);
			num_allele_2.push_back(0);
			freq_allele_1.push_back(0);
			freq_allele_2.push_back(0);
		}

		int maternal_1, maternal_2, paternal_1, paternal_2;

		for (fff = 0; fff < pop_size; fff++)
		{
			maternal_1 = adults[fff].maternal[chromosome_1].loci[allele_1];
			maternal_2 = adults[fff].maternal[chromosome_2].loci[allele_2];
			paternal_1 = adults[fff].paternal[chromosome_1].loci[allele_1];
			paternal_2 = adults[fff].paternal[chromosome_2].loci[allele_2];

			num_allele_1[maternal_1]++;
			num_allele_1[paternal_1]++;
			num_allele_2[maternal_2]++;
			num_allele_2[paternal_2]++;
			joint_alleles[maternal_1][maternal_2]++;
			joint_alleles[paternal_1][paternal_2]++;
		}//end for fff < PopulationSize

		int major_allele_1 = 0;
		int major_allele_2 = 0;
		for (al = 0; al < num_alleles; al++)
		{
			freq_allele_1[al] = (num_allele_1[al]) / (2 * dNadult);
			freq_allele_2[al] = (num_allele_2[al]) / (2 * dNadult);
			if (freq_allele_1[al] > freq_allele_1[major_allele_1])
				major_allele_1 = al;
			if (freq_allele_2[al] > freq_allele_2[major_allele_2])
				major_allele_2 = al;
			for (f = 0; f < num_alleles; f++)
				joint_alleles[al][f] = joint_alleles[al][f] / (2 * dNadult);
		}
		primary_frequency[chromosome_1].locus_emulator[allele_1] = freq_allele_1[major_allele_1];
		major_allele[chromosome_1].locus_emulator[allele_1] = major_allele_1;
		primary_frequency[chromosome_2].locus_emulator[allele_2] = freq_allele_2[major_allele_2];
		major_allele[chromosome_2].locus_emulator[allele_2] = major_allele_2;


		if (freq_allele_1[major_allele_1] != 1 && freq_allele_2[major_allele_1] != 1)//if it's polymorphic
		{
			d_allele_avgs = 0;
			count = 0;
			for (f = 0; f < num_alleles; f++)
			{
				for (ff = 0; ff < num_alleles; ff++)
				{
					if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
					{
						D[f][ff] = joint_alleles[f][ff] - freq_allele_1[f] * freq_allele_2[ff];
						d_allele_avgs = d_allele_avgs + fabs(D[f][ff]);
						count++;
						if (D[f][ff] < 0)
							Dmax[f][ff] = min(freq_allele_1[f] * freq_allele_2[ff], (1 - freq_allele_1[f])*(1 - freq_allele_2[ff]));
						else
							Dmax[f][ff] = min((1 - freq_allele_1[f])*freq_allele_2[ff], freq_allele_1[f] * (1 - freq_allele_2[ff]));
					}
				}
			}
			double dcount = count;
			d_allele_avgs = d_allele_avgs / dcount;

			Dprime = 0;
			bool decentDmax = true;
			for (f = 0; f < num_alleles; f++)
			{
				for (ff = 0; ff < num_alleles; ff++)
				{
					if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
					{
						if (Dmax[f][ff] > 0)
							Dprime = Dprime + freq_allele_1[f] * freq_allele_2[ff] * fabs(D[f][ff]) / Dmax[f][ff];
						else
							decentDmax = false;

					}
				}
			}
			if (!decentDmax)
				Dprime = -5;

		}

		result.d = d_allele_avgs;
		result.dprime = Dprime;



		for (f = 0; f < num_alleles; f++)
			delete[] joint_alleles[f];
		for (f = 0; f < num_alleles; f++)
			delete[] D[f];
		for (f = 0; f < num_alleles; f++)
			delete[] Dmax[f];
		delete[] joint_alleles;
		delete[] D;
		delete[] Dmax;

		return result;
	}//end Adult Pop LD

	void standardize_genotypes()
	{
		int j, jj, jjj;
		vector<double> alleles;
		double mean, sd, constant;
		int num_alleles = 0;
		mean = sd = 0;
		for (j = 0; j < pop_size; j++)
			mean = mean + adults[j].genotype;
		mean = mean / pop_size;
		for (j = 0; j < pop_size; j++)
			sd = sd + (adults[j].genotype - mean)*(adults[j].genotype - mean);
		sd = sqrt(sd / pop_size);
		constant = mean / (num_qtl * 2 * num_chrom);
		for (j = 0; j < pop_size; j++)
		{
			adults[j].phenotype = adults[j].phenotype - adults[j].genotype;
			adults[j].genotype = 0;
			for (jj = 0; jj < num_chrom; jj++)
			{
				for (jjj = 0; jjj < num_qtl; jjj++)
				{
					adults[j].paternal[jj].allelic_effects[jjj] = (adults[j].paternal[jj].allelic_effects[jjj] - constant) / sd;
					adults[j].maternal[jj].allelic_effects[jjj] = (adults[j].maternal[jj].allelic_effects[jjj] - constant) / sd;
					adults[j].genotype = adults[j].genotype + adults[j].paternal[jj].allelic_effects[jjj];
				}
			}
			adults[j].phenotype = adults[j].genotype + adults[j].phenotype;
		}
		mean = sd = 0;
		for (j = 0; j < pop_size; j++)
			mean = mean + adults[j].genotype;
		mean = mean / pop_size;
		for (j = 0; j < pop_size; j++)
			sd = sd + (adults[j].genotype - mean)*(adults[j].genotype - mean);
		sd = sqrt(sd / pop_size);
	}

	void dads_genotype(int adult_index, int off_index, int dad_index)
	{
		int j, jj, jjj;

		for (j = 0; j < num_chrom; j++)
		{
			for (jj = 0; jj < num_markers; jj++)
			{
				inferred_dads[dad_index][j].loci[jj] = offspring[off_index].paternal[j].loci[jj];
			}
			for (jj = 0; jj < num_qtl; jj++)
			{
				inferred_dads[dad_index][j].allelic_effects[jj] = offspring[off_index].paternal[j].allelic_effects[jj];
			}
		}
	}
	
	void infer_genotype(int adult_index, int off_index, int dad_index)
	{
		int j, jj, jjj;
		for (j = 0; j < num_chrom; j++)
		{
			for (jj = 0; jj < num_markers; jj++)
			{//sample the way I sampled the actual population
				vector <int> mom_allele;
				if (adults[adult_index].paternal[j].loci[jj] == offspring[off_index].paternal[j].loci[jj])
					mom_allele.push_back(offspring[off_index].paternal[j].loci[jj]);
				if (adults[adult_index].paternal[j].loci[jj] == offspring[off_index].maternal[j].loci[jj])
					mom_allele.push_back(offspring[off_index].maternal[j].loci[jj]);
				if (adults[adult_index].maternal[j].loci[jj] == offspring[off_index].paternal[j].loci[jj])
					mom_allele.push_back(offspring[off_index].paternal[j].loci[jj]);
				if (adults[adult_index].maternal[j].loci[jj] == offspring[off_index].maternal[j].loci[jj])
					mom_allele.push_back(offspring[off_index].maternal[j].loci[jj]);

				if (mom_allele.size() > 1)
				{
					vector<int>::iterator it;
					it = unique(mom_allele.begin(), mom_allele.end());
					mom_allele.resize(distance(mom_allele.begin(), it));
				}

				if (mom_allele.size() == 1)//otherwise it's not an informative locus
				{
					inferred_dads[dad_index][j].loci[jj] = mom_allele[0];
				}
			}
			for (jj = 0; jj < num_qtl; jj++)
			{
				inferred_dads[dad_index][j].allelic_effects[jj] = offspring[off_index].paternal[j].allelic_effects[jj];
			}
		}
	}

	void calculate_pop_numbers()
	{
		int j;
		pop_size = num_fem = num_mal = num_off = 0;
		for (j = 0; j < adults.size(); j++)
		{
			if (adults[j].alive)
			{
				if (adults[j].female)
					num_fem++;
				else
					num_mal++;
				pop_size++;
			}
		}
		for (j = 0; j < offspring.size(); j++)
		{
			if (offspring[j].alive)
				num_off++;
		}
	}

	void sample_pop()
	{
		int j, jj, jjj, off_counter, fem_counter, mal_counter;
		vector<bool> taken;
		int rand1, rand2;
		
		calculate_pop_numbers();

		if (adult_samplesize > pop_size)
			adult_samplesize = pop_size;
		if (female_samplesize > num_fem)
			female_samplesize = num_fem;
		if (male_samplesize > num_mal)
			male_samplesize = num_mal;
		if (offspring_samplesize > num_off)
			offspring_samplesize = num_off;
		dad_samplesize = offspring_samplesize;
		for (j = 0; j < pop_size; j++)
			taken.push_back(false);

		//reset the vectors
		sampled_adults.resize(0);
		sampled_females.resize(0);
		sampled_males.resize(0);
		sampled_dads.resize(0);
		sampled_off.resize(0);
		//from a single collection of adults, we'll fill in all of the other things
		off_counter = fem_counter = mal_counter = 0;
		for (j = 0; j < adult_samplesize; j++)
		{
			rand1 = randnum(pop_size);
			if (taken[rand1] == false)
			{
				if (adults[rand1].alive)//only sample if it's alive
				{
					sampled_adults.push_back(rand1);
					if (adults[sampled_adults[j]].female)
					{
						if (fem_counter < female_samplesize)
						{
							sampled_females.push_back(rand1);
							if (off_counter < offspring_samplesize)
							{
								rand2 = randnum(adults[sampled_females[fem_counter]].offspring_index.size());
								sampled_off.push_back(adults[sampled_females[fem_counter]].offspring_index[rand2]);
								//infer dad's genotype
								infer_genotype(sampled_adults[j], sampled_off[off_counter], off_counter);
								sampled_dads.push_back(off_counter);
							}
						}
						off_counter++;
						fem_counter++;
					}
					else
					{
						if (mal_counter < male_samplesize)
							sampled_males.push_back(rand1);
					}
				}
			}
			else
			{
				while (taken[rand1] == true)
					rand1 = randnum(pop_size);
				if (taken[rand1] == false)
				{
					if (adults[sampled_adults[j]].alive)//only sample live individuals
					{
						if (adults[rand1].alive)//only sample if it's alive
						{
							sampled_adults.push_back(rand1);
							if (adults[sampled_adults[j]].female)
							{
								if (fem_counter < female_samplesize)
								{
									sampled_females.push_back(rand1);
									if (off_counter < offspring_samplesize)
									{
										rand2 = randnum(adults[sampled_females[fem_counter]].offspring_index.size());
										sampled_off.push_back(adults[sampled_females[fem_counter]].offspring_index[rand2]);
										//infer dad's genotype
										infer_genotype(sampled_adults[j], sampled_off[off_counter], off_counter);
									}
								}
								off_counter++;
								fem_counter++;
							}
							else
							{
								if (mal_counter < male_samplesize)
									sampled_males.push_back(rand1);
							}
						}
					}
				}//if you find one
			}//else
		}
	}//sample pop

	void CalculateAdultMarkerAlleleFreqs(int whichchromosome, int marker, int ref, locus_statistics &fsts, vector<int>&sampled, bool pop1, bool includeDead)
	{
		int caa;
		vector<double> tempallelefreq;
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq.push_back(0);

		double allelecounter = 0;
		for (caa = 0; caa < sampled.size(); caa++)
		{
			if (adults[sampled[caa]].alive || includeDead)
			{
				tempallelefreq[adults[sampled[caa]].maternal[whichchromosome].loci[marker]]++;
				tempallelefreq[adults[sampled[caa]].paternal[whichchromosome].loci[marker]]++;
				allelecounter = allelecounter + 2;
			}
		}
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq[caa] = tempallelefreq[caa] / allelecounter;
		int allelesperlocus = 0;
		for (caa = 0; caa < num_alleles; caa++)
		{
			if (pop1)
				fsts.pop1_freqs[whichchromosome].freqs[ref][caa] = tempallelefreq[caa];
			else
				fsts.pop2_freqs[whichchromosome].freqs[ref][caa] = tempallelefreq[caa];
			if (tempallelefreq[caa] != 0)
				allelesperlocus++;
		}
		fsts.num_alleles[whichchromosome][ref] = allelesperlocus;
		
	}//end CalculateAdultMarkerAlleleFreqs

	void CalculateProgenyMarkerAlleleFreqs(int whichchromosome, int marker, int ref, bool includeDead)
	{
		int caa;
		vector<double> tempallelefreq;
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq.push_back(0);

		double allelecounter = 0;
		for (caa = 0; caa < offspring_samplesize; caa++)
		{
			if (offspring[sampled_off[caa]].alive || includeDead)
			{
				tempallelefreq[offspring[sampled_off[caa]].maternal[whichchromosome].loci[marker]]++;
				tempallelefreq[offspring[sampled_off[caa]].paternal[whichchromosome].loci[marker]]++;
				allelecounter = allelecounter + 2;
			}
		}
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq[caa] = tempallelefreq[caa] / allelecounter;
		int allelesperlocus = 0;
		for (caa = 0; caa < num_alleles; caa++)
		{
			adult_offspring.pop2_freqs[whichchromosome].freqs[ref][caa] = tempallelefreq[caa];
			if (tempallelefreq[caa] != 0)
				allelesperlocus++;
		}
		adult_offspring.num_alleles[whichchromosome][ref] = allelesperlocus;
	}

	void CalculateDadsMarkerAlleleFreqs(int whichchromosome, int marker, int ref)
	{
		int caa;
		vector<double> tempallelefreq;
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq.push_back(0);
		double allelecounter = 0;
		for (caa = 0; caa < inferred_dads.size(); caa++)
		{
			tempallelefreq[inferred_dads[caa][whichchromosome].loci[marker]]++;
			allelecounter = allelecounter++;
		}
		for (caa = 0; caa < num_alleles; caa++)
			tempallelefreq[caa] = tempallelefreq[caa] / allelecounter;
		int allelesperlocus = 0;
		for (caa = 0; caa < num_alleles; caa++)
		{
			males_dads.pop2_freqs[whichchromosome].freqs[ref][caa] = tempallelefreq[caa];
			if (tempallelefreq[caa] != 0)
				allelesperlocus++;
		}
		males_dads.num_alleles[whichchromosome][ref] = allelesperlocus;
	}

	double determine_major_af(vector<double> allele_freqs)
	{
		int j;
		double max = 0;
		for (j = 0; j < num_alleles; j++)
		{
			if (allele_freqs[j] > max)
				max = allele_freqs[j];
		}
		return max;
	}

	void calculate_fsts(locus_statistics& fsts, vector<int> pop1, vector<int> pop2, bool pop2_off, bool pop2_dads, string file_name)
	{
		ofstream fst_out;
		string fst_out_name = file_name + ".txt";
		fst_out.open(fst_out_name);
		fst_out << "Chrom\tLocus\tPop1N\tPop2N\tPop1AF\tPop1AlleleN\tPop2AF\tPop2AlleleN\tPop1Hs\tPop2Hs\tAvgHs\tHt\tFst";
		vector<double> averageallelefreqs;
		for (int j = 0; j < num_alleles; j++)
			averageallelefreqs.push_back(0);
		int locus;
		double hs1, hs2, ht, fst;
		int whichallele, whichchromosome, marker;
		for (whichchromosome = 0; whichchromosome < num_chrom; whichchromosome++)
		{
			locus = 0;
			for (marker = 0; marker < num_markers; marker++)
			{
				fsts.locus_id[whichchromosome].locus_emulator[marker] = marker;
				CalculateAdultMarkerAlleleFreqs(whichchromosome, marker, locus, fsts, pop1, true, false);
				fst_out << '\n' << whichchromosome << '\t' << marker << '\t' << pop1.size() << '\t' << pop2.size()
					<< '\t' << determine_major_af(fsts.pop1_freqs[whichchromosome].freqs[locus])
					<< '\t' << fsts.num_alleles[whichchromosome][marker];
				if (pop2_off)
					CalculateProgenyMarkerAlleleFreqs(whichchromosome, marker, locus, false);
				if (pop2_dads)
					CalculateDadsMarkerAlleleFreqs(whichchromosome, marker, locus);
				if (pop2_off == false && pop2_dads == false)
					CalculateAdultMarkerAlleleFreqs(whichchromosome, marker, locus, fsts, pop2, false, false);
				for (whichallele = 0; whichallele < num_alleles; whichallele++)
				{
					averageallelefreqs[whichallele] = (fsts.pop1_freqs[whichchromosome].freqs[locus][whichallele] +
						fsts.pop2_freqs[whichchromosome].freqs[locus][whichallele]) / 2;
				}
				fst_out << '\t' << determine_major_af(fsts.pop2_freqs[whichchromosome].freqs[locus])
					<< '\t' << fsts.num_alleles[whichchromosome][locus];
				hs1 = 1;
				hs2 = 1;
				ht = 1;

				for (whichallele = 0; whichallele < num_alleles; whichallele++)
				{
					hs1 = hs1 - (fsts.pop1_freqs[whichchromosome].freqs[locus][whichallele] * fsts.pop1_freqs[whichchromosome].freqs[locus][whichallele]);
					hs2 = hs2 - (fsts.pop2_freqs[whichchromosome].freqs[locus][whichallele] * fsts.pop2_freqs[whichchromosome].freqs[locus][whichallele]);
					ht = ht - (averageallelefreqs[whichallele] * averageallelefreqs[whichallele]);
				}
				if (ht > 0)
					fst = 1 - (hs1 + hs2) / (2 * ht);
				else
					fst = -1.0;
				fst_out << '\t' << hs1 << '\t' << hs2 << '\t' << (hs1 + hs2) / 2 << '\t' << ht << '\t' << fst;
				fsts.fst[whichchromosome].locus_emulator[locus] = fst;
				fsts.exp_het[whichchromosome].locus_emulator[locus] = ht;
				locus++;
			}
		}
		fst_out.close();
	}//end calculate_fsts

	
};//end Population