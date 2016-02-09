//Author: Sarah P. Flanagan (sflanagan@bio.tamu.edu)
//Date: 4 December 2015
//Purpose: Generate a null distribution for a selection components analysis
//SCA compares adults to offspring, males to females, and inferred parents to the population

#include "population.h"	
#include "random_numbers.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <array>

using namespace std;

int main()
{
	int num_reps = 1;
	string base_name = "knowndist.ss0.2alleles";
	bool known_qtl = false;
	bool empirical_afs = true;
	int end, generations, reps, i, ii, iii, ld_count;
	population pop;
	ld_info returned_data;
	double mean_dp_ldistchrom;
	string ld_out_name = base_name + ".ld_out.txt";
	string fst_out_name = base_name + ".fst_out.txt";
	ofstream fst_out;
	ofstream ld_out;
	
	ld_out.open(ld_out_name);
	ld_out << "Reps\tGens\tMeanD'LD\tMeanPairwiseD\tMeanLD";
	for (reps = 0; reps < num_reps; reps++)
	{
		pop.set_parameters();
		pop.initialize(known_qtl, empirical_afs);
		for (generations = 0; generations < pop.num_gen; generations++)
		{
			pop.determine_pop_size();
			pop.mating(0);
			if (pop.extinct == true)
			{
				cout << "\nExtinct generation " << generations + 1 << "\n";
				break;
			}
			pop.mutation();
			pop.viability_selection(0);
			pop.density_regulation();
			pop.calc_mean_trait_values();
			//ld summary stats
			double mean_new_longdist_d, mean_new_longdist_dprime;
			mean_new_longdist_d = mean_new_longdist_dprime = 0;
			ld_count = 0;
			while (ld_count < pop.num_ld)
			{
				returned_data = pop.adult_pop_ld(randnum(pop.num_chrom), randnum(pop.num_chrom), randnum(pop.num_markers), randnum(pop.num_markers));
				if (returned_data.dprime != -5)
				{
					mean_new_longdist_d = mean_new_longdist_d + returned_data.d;
					mean_new_longdist_dprime = mean_new_longdist_dprime + returned_data.dprime;
					ld_count++;
				}
			}
			mean_new_longdist_d = mean_new_longdist_d / pop.num_ld;
			mean_new_longdist_dprime = mean_new_longdist_dprime / pop.num_ld;
			for (int c = 0; c < pop.num_chrom; c++)
			{
				for (int cc = 0; cc < pop.num_markers - 1; cc++)
				{//Pairwise
					returned_data = pop.adult_pop_ld(c, c, cc, cc + 1);
					pop.avg_pairwise_d = pop.avg_pairwise_d + returned_data.d;
					pop.avg_ld = pop.avg_ld + returned_data.dprime;
				}
				//adams_avg_dprime[c] = 0;
				pop.pop_ld[c].d = 0;
				ld_count = 0;
				pop.pop_ld[c].dprime = 0;
				mean_dp_ldistchrom = 0;
				while (ld_count < 100)
				{//random within-chromosomes
					returned_data = pop.adult_pop_ld(c, c, randnum(pop.num_markers), randnum(pop.num_markers));
					if (returned_data.dprime != -5)
					{
						pop.pop_ld[c].d = pop.pop_ld[c].d + returned_data.d;
						pop.pop_ld[c].dprime = pop.pop_ld[c].dprime + returned_data.dprime;
						mean_dp_ldistchrom = mean_dp_ldistchrom + returned_data.dprime;
						ld_count++;
					}
				}
				pop.pop_ld[c].d = pop.pop_ld[c].d / 100;
				pop.pop_ld[c].dprime = pop.pop_ld[c].dprime / 100;
			}
			mean_dp_ldistchrom = mean_dp_ldistchrom / (pop.num_chrom * 100);
			pop.avg_ld = pop.avg_ld / (pop.num_chrom*pop.num_markers);
			pop.avg_pairwise_d = pop.avg_pairwise_d / (pop.num_chrom*pop.num_markers);
			ld_out << '\n' << reps << '\t' << generations << '\t' << mean_dp_ldistchrom << '\t' << pop.avg_pairwise_d << '\t' << pop.avg_ld;
			if ((generations + 1) % 100 == 0)
				cout << "Generation " << generations + 1 << " complete.\t";
		}
	}
	cout << '\n';

	if (!pop.extinct)
	{
		pop.pop_size = pop.determine_pop_size();
		if (pop.pop_size == 0)
			pop.extinct = true;
		else
		{
			cout << "\nSampling the population.\n";
			pop.standardize_genotypes();
			pop.mating(0);
			pop.mutation();
			pop.sample_pop();
			ofstream sampled_inds;
			sampled_inds.open("sampled.inds.txt");
			sampled_inds << "Adults";
			for (int sampled = 0; sampled < pop.sampled_adults.size(); sampled++)
				sampled_inds << '\t' << pop.sampled_adults[sampled];
			sampled_inds << "\nOffspring";
			for (int sampled = 0; sampled < pop.sampled_off.size(); sampled++)
				sampled_inds << '\t' << pop.sampled_off[sampled];
			sampled_inds << "\nMales";
			for (int sampled = 0; sampled < pop.sampled_males.size(); sampled++)
				sampled_inds << '\t' << pop.sampled_males[sampled];
			sampled_inds << "\nFemales";
			for (int sampled = 0; sampled < pop.sampled_females.size(); sampled++)
				sampled_inds << '\t' << pop.sampled_females[sampled];
			sampled_inds << "\nDads";
			for (int sampled = 0; sampled < pop.sampled_dads.size(); sampled++)
				sampled_inds << '\t' << pop.sampled_dads[sampled];
			sampled_inds.close();

			pop.calculate_fsts(pop.adult_offspring, pop.sampled_adults, pop.sampled_off, true, false, base_name + "ao");
			pop.calculate_fsts(pop.male_female, pop.sampled_males, pop.sampled_females, false, false, base_name + "mf");
			pop.calculate_fsts(pop.males_dads, pop.sampled_males, pop.sampled_dads, false, true, base_name + "gp");
			//Calculate LD
			for (int c = 0; c < pop.num_chrom; c++)
			{
				for (int cc = 0; cc < pop.num_markers - 1; cc++)
				{//Pairwise
					returned_data = pop.adult_pop_ld(c, c, cc, cc + 1);
					pop.avg_pairwise_d = pop.avg_pairwise_d + returned_data.d;
					pop.avg_ld = pop.avg_ld + returned_data.dprime;
				}
			}
			pop.avg_ld = pop.avg_ld / (pop.num_chrom*pop.num_markers);
			pop.avg_pairwise_d = pop.avg_pairwise_d / (pop.num_chrom*pop.num_markers);
			fst_out.open(fst_out_name);
			fst_out << "Locus\tAOFst\tMFFst\tMDFst\tAdultAF\tAdultN\tOffAF\tOffN\tMaleAF\tMaleN\tFemAF\tFemN\tDadAF\tDadN";
			for (i = 0; i < pop.num_chrom; i++)
			{
				for (ii = 0; ii < pop.num_markers; ii++)
				{
					fst_out << '\n' << pop.adult_offspring.locus_id[i].locus_emulator[ii];
					fst_out << '\t' << pop.adult_offspring.fst[i].locus_emulator[ii];
					fst_out << '\t' << pop.male_female.fst[i].locus_emulator[ii];
					fst_out << '\t' << pop.males_dads.fst[i].locus_emulator[ii];
					fst_out << '\t' << pop.determine_major_af(pop.adult_offspring.pop1_freqs[i].freqs[ii]);
					fst_out << '\t' << pop.adult_offspring.num_alleles[i][ii];
					fst_out << '\t' << pop.determine_major_af(pop.adult_offspring.pop2_freqs[i].freqs[ii]);
					fst_out << '\t' << pop.adult_offspring.num_alleles[i][ii];
					fst_out << '\t' << pop.determine_major_af(pop.male_female.pop1_freqs[i].freqs[ii]);
					fst_out << '\t' << pop.male_female.num_alleles[i][ii];
					fst_out << '\t' << pop.determine_major_af(pop.male_female.pop2_freqs[i].freqs[ii]); 
					fst_out << '\t' << pop.male_female.num_alleles[i][ii];
					fst_out << '\t' << pop.determine_major_af(pop.males_dads.pop2_freqs[i].freqs[ii]); 
					fst_out << '\t' << pop.males_dads.num_alleles[i][ii];
				}
			}
			fst_out.close();
		}
	}
	ld_out.close();

	cout << "\nDone! Input Integer to Quit.\n";
	cin >> end;
	return 0;
}