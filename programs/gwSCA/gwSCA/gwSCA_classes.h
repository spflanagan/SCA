#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "math.h"

using namespace std;

void FileTest(ifstream& file, string filename)
{
	cout << filename;
	if (file.is_open())
		cout << " open\n";
	else
	{
		while (!file.is_open())
		{
			cout << " not open. Please re-enter filename: ";
			getline(cin, filename, '\n');
			file.open(filename);
		}
	}

}

istream& universal_getline(istream& is, string& t)
{
	//this code is adapted from a post on stackoverflow:
	// http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	//written by user763305
	t.clear();
	istream::sentry se(is, true);
	streambuf* sb = is.rdbuf();//sets pointer to stream buffer object

	for (;;)
	{
		int c = sb->sbumpc();//get current character and advance to the next position
		switch (c)//tests for equality against a list of variables (like multiple if statements)
		{
		case '\n'://if the next character is '\n', return the line
			return is;
		case '\r'://if the character is '\n', see if the next one is '\n' or return the line
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:// Also handle the case when the last line has no line ending
			if (t.empty())//if there's nothing there, set it to be the end of file
				is.setstate(ios::eofbit);//set it to be the end of the file and return it
			return is;
		default://if none of the above, continue on.
			t += (char)c;
		}
	}

}

string find_and_replace(string &s, string toReplace, string replaceWith)
{
	size_t pos = 0;
	while ((pos = s.find(toReplace, pos)) != string::npos) {
		s.replace(pos, toReplace.length(), replaceWith);
		pos += replaceWith.length();
	}
	return s;
}

class Cutoffs
{
public:
	double critical_99;
	double critical_98;
	double critical_95;
	double critical_90;
	double critical_80;
};
class Locus_info
{
public:
	string snp_name;
	vector <char> snps;
	vector <double> snp_freq;
	double distance;

	Locus_info()
	{
		snp_name = "NULL";
		snps.push_back('n');
		snp_freq.push_back(0);
	}
};

class Chromosome
{
public:
	string ID;
	vector <Locus_info> loci;

	Chromosome()
	{
		loci.push_back(Locus_info());
		ID = "NULL";
	}

public:
	void update_chromosome(char input, int locus)
	{
		bool found = false;
		size_t i, ii;

		if (input != '0')
		{
			for (i = 0; i < loci[locus].snps.size(); i++)
			{
				if (input == loci[locus].snps[i])
				{
					found = true;
				}
			}
			if (found == false)
			{
				loci[locus].snps.push_back(input);
				loci[locus].snp_freq.push_back(0);
			}
		}
	}
	void remove_defaults()
	{
		size_t ii, iii;
		for (ii = 0; ii < loci.size(); ii++)
		{
			for (iii = 0; iii < loci[ii].snps.size(); iii++)
			{
				if (loci[ii].snps[iii] == 'n')
				{
					int this_locus = iii;
					loci[ii].snps.erase(loci[ii].snps.begin() + this_locus);
					loci[ii].snp_freq.erase(loci[ii].snp_freq.begin() + this_locus);

				}
			}
			if (loci[ii].snps.size() == 0)
				loci.erase(loci.begin() + ii);
		}
	}
};

class Genome
{
public:
	vector <Chromosome> chrom_set;

	Genome(const vector<Chromosome> &ref)
	{
		size_t l, ll, lll;
		for (l = 0; l < ref.size(); l++)
		{
			chrom_set.push_back(Chromosome());
			for (ll = 0; ll < ref[l].loci.size(); ll++)
			{
				chrom_set[l].loci.push_back(Locus_info());
				chrom_set[l].loci[ll].snp_name = ref[l].loci[ll].snp_name;
				for (lll = 0; lll < ref[l].loci[ll].snp_freq.size(); lll++)
				{
					chrom_set[l].loci[ll].snp_freq.push_back(0);
					chrom_set[l].loci[ll].snps.push_back(ref[l].loci[ll].snps[lll]);
				}
			}

		}
	}
};

class Individual
{
public:
	string ID;
	bool female;
	vector <int> phenotype1;
	vector <int> phenotype2;
	vector <int> year;
	vector <bool> phenotype3;
	vector <Chromosome> allele1;
	vector <Chromosome> allele2;

	Individual()
	{
		ID = "NULL";
		female = true;
		phenotype1.push_back(0);
		phenotype2.push_back(0);
		year.push_back(0);
		phenotype3.push_back(0);
		allele1.push_back(Chromosome());
		allele2.push_back(Chromosome());
	}
};//end Individual

class Population
{
public:
	string name;
	vector <Individual> adults;
	vector <int> all_inds;
	vector < vector <int> > subsets;

	Population()
	{
		name = "";
		adults.push_back(Individual());
	}

public:
	void remove_defaults()
	{
		size_t ii, iii;
		for (ii = 0; ii < adults.size(); ii++)
		{
			if (adults[ii].ID == "NULL")
			{
				int this_locus = ii;
				adults.erase(adults.begin() + this_locus);
			}
		}
	}

	void establish_subsets(int phen_no, bool heirarchical, int thresh)
	{
		int count, last, subset_no;
		size_t s, ss, sss;
		for (s = 0; s < adults.size(); s++)
		{
			all_inds.push_back(s);
		}
		if (phen_no == 1)
		{
			vector <int> temp_subsets;
			bool found;
			temp_subsets.push_back(adults[0].phenotype1[0]);
			count = 0;
			for (s = 0; s < adults.size(); s++)
			{
				for (ss = 0; ss < adults[s].phenotype1.size(); ss++)
				{
					if (adults[s].phenotype1[ss] != temp_subsets[count])
					{
						found = false;
						for (sss = 0; sss < temp_subsets.size(); sss++)
						{
							if (temp_subsets[sss] == adults[s].phenotype1[ss])
								found = true;
						}
						if (found == false)
						{
							count++;
							temp_subsets.push_back(adults[s].phenotype1[ss]);
						}
					}
				}
			}
			subset_no = temp_subsets.size();
			//	cout << "number of subsets: " << subset_no << '\n';
			for (s = 0; s < subset_no; s++)
			{
				vector<int> this_subset;
				subsets.push_back(this_subset);
			}
			if (heirarchical == true)//then the maximum value trumps lower ones. (3>2>1>0)
			{
				int last_value;
				for (s = 0; s < adults.size(); s++)
				{
					last_value = 0;
					for (ss = 0; ss < adults[s].phenotype1.size(); ss++)
					{
						if (adults[s].phenotype1[ss] > last_value)
							last_value = adults[s].phenotype1[ss];
					}
					if (last_value < subsets.size())
						subsets.at(last_value).push_back(s);
					else
						cout << "Adult " << s << " has a group value of " << last_value << ", which is larger than " << subset_no
						<< ", the number of groups.\n";
				}
			}
			else//it's the first value
			{
				for (s = 0; s < adults.size(); s++)
				{
					all_inds.push_back(s);
					if (adults[s].phenotype1[0] < subsets.size())
						subsets.at(adults[s].phenotype1[ss]).push_back(s);
					else
						cout << "Adult " << s << " has a group value of " << adults[s].phenotype1[0] << ", which is larger than " << subset_no
						<< ", the number of groups.\n";
					/*for (ss = 0; ss < adults[ss].phenotype1.size(); ss++)
					{
					if (adults[s].phenotype1[ss] < subsets.size())
					subsets.at(adults[s].phenotype1[ss]).push_back(s);
					else
					cout << "Adult " << s << " has a group value of " << adults[s].phenotype1[ss] << ", which is larger than " << subset_no
					<< ", the number of groups.\n";
					}*/
				}
			}
			for (s = 0; s < subset_no; s++)
			{
				cout << name << " has " << subsets[s].size() << " individuals of group " << s << "\n";
			}
		}//phen1
		if (phen_no == 2)//use the threshold
		{
			subset_no = 2; //sorting into those below threshold and those above
			bool below;
			for (s = 0; s < subset_no; s++)
			{
				vector<int> this_subset;
				subsets.push_back(this_subset);
			}
			for (s = 0; s < adults.size(); s++)
			{
				below = true;
				for (ss = 0; ss < adults[s].phenotype2.size(); ss++)
				{
					if (adults[s].phenotype2[ss] > thresh)
					{
						below = false;
					}
				}
				if (below)
				{
					subsets[0].push_back(s);
				}
				else
				{
					subsets[1].push_back(s);
				}
			}
			for (s = 0; s < subset_no; s++)
			{
				cout << name << " has " << subsets[s].size() << " individuals of group " << s << "\n";
			}
		}//end phen2
		if (phen_no == 3)//it's boolean
		{
			subset_no = 2;
			bool zero;
			for (s = 0; s < subset_no; s++)
			{
				vector<int> this_subset;
				subsets.push_back(this_subset);
			}
			for (s = 0; s < adults.size(); s++)
			{
				zero = true;
				for (ss = 0; ss < adults[s].phenotype3.size(); ss++)
				{
					if (adults[s].phenotype3[ss])
					{
						zero = false;
					}
				}
				if (zero)
				{
					subsets[0].push_back(s);
				}
				if (!zero)
				{
					subsets[1].push_back(s);
				}
			}
			for (s = 0; s < subset_no; s++)
			{
				cout << name << " has " << subsets[s].size() << " individuals of group " << s << "\n";
			}
		}//end phen3
		if (phen_no == 4)//use 1 and 3
		{
			vector <int> temp_subsets;
			bool found;
			int max_phen1 = 0;
			temp_subsets.push_back(adults[0].phenotype1[0]);
			count = 0;
			for (s = 0; s < adults.size(); s++)//must be hierarchical
			{
				for (ss = 0; ss < adults[s].phenotype1.size(); ss++)
				{
					if (adults[s].phenotype1[ss] > max_phen1)
					{
						max_phen1 = adults[s].phenotype1[ss];
					}
				}
			}
			subset_no = 2;
			for (s = 0; s < 2; s++)
			{
				vector<int> this_subset;
				subsets.push_back(this_subset);
			}

			for (s = 0; s < adults.size(); s++)
			{
				found = false;
				for (ss = 0; ss < adults[s].year.size() - 1; ss++)
				{
					if (adults[s].phenotype3[ss] == 1)
					{
						if (adults[s].phenotype1[ss + 1] == max_phen1)
							found = true;
					}
				}
				if (found)
				{
					subsets[1].push_back(s);
				}
				else
				{
					subsets[0].push_back(s);
				}
			}
			for (s = 0; s < subset_no; s++)
			{
				cout << name << " has " << subsets[s].size() << " individuals of group " << s << "\n";
			}
		}//end phen 4
	}//end establish_subsets

	Genome calc_allele_freq(vector <int> &subset, vector <Chromosome> &ref)
	{

		//cout << "calculating allele frequencies\n";
		size_t i, ii, iii, s_chr;
		Genome allele_freqs(ref);
		double num_adults = subset.size();
		int count = 0;
		double avg_freq = 0;
		double avg_count = 0;

		for (s_chr = 0; s_chr < ref.size(); s_chr++)
		{
			for (ii = 0; ii < ref[s_chr].loci.size(); ii++)
			{
				for (iii = 0; iii < ref[s_chr].loci[ii].snp_freq.size(); iii++)
				{
					allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] = 0;
					allele_freqs.chrom_set[s_chr].loci[ii].snps[iii] = ref[s_chr].loci[ii].snps[iii];
				}
			}
		}
		for (i = 0; i < subset.size(); i++)//loop through all of the selected adults
		{
			for (s_chr = 0; s_chr < ref.size(); s_chr++)
			{
				for (ii = 0; ii < ref[s_chr].loci.size(); ii++)
				{
					for (iii = 0; iii < ref[s_chr].loci[ii].snps.size(); iii++)
					{
						if (adults[subset[i]].allele1[s_chr].loci[ii].snps[0] == ref[s_chr].loci[ii].snps[iii])
						{
							allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] = allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] + 1;
							count++;
						}

						if (adults[subset[i]].allele2[s_chr].loci[ii].snps[0] == ref[s_chr].loci[ii].snps[iii])
						{
							allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] = allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] + 1;
							count++;
						}
					}
				}
			}
		}
		for (s_chr = 0; s_chr < ref.size(); s_chr++)
		{
			for (ii = 0; ii < ref[s_chr].loci.size(); ii++)
			{
				for (iii = 0; iii < ref[s_chr].loci[ii].snps.size(); iii++)
				{
					//	cout << allele_freqs.chrom_set[s_chr].loci[ii].snps[iii] << '\t' <<  allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] << '\t';
					allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] = allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] / (num_adults * 2);
					//	cout << allele_freqs.chrom_set[s_chr].loci[ii].snp_freq[iii] << '\n';
				}
			}
		}
		return allele_freqs;
	}//end calc_allele_freq

	Genome calc_obs_het(vector<int> &subset, vector <Chromosome> &ref)
	{
		Genome obs_het(ref);
		size_t i, ii, iii, s_chr;
		double num_adults = subset.size();

		for (i = 0; i < subset.size(); i++)//loop through all of the selected adults
		{
			for (s_chr = 0; s_chr < ref.size(); s_chr++)
			{
				for (ii = 0; ii < ref[s_chr].loci.size(); ii++)
				{
					if (adults[subset[i]].allele1[s_chr].loci[ii].snps[0] == adults[subset[i]].allele2[s_chr].loci[ii].snps[0])
					{
						//cout << "het found!\t";
						for (iii = 0; iii < ref[s_chr].loci[ii].snps.size(); iii++)
						{
							if (adults[subset[i]].allele1[s_chr].loci[ii].snps[0] == ref[s_chr].loci[ii].snps[iii])
								obs_het.chrom_set[s_chr].loci[ii].snp_freq[iii] = obs_het.chrom_set[s_chr].loci[ii].snp_freq[iii] + 1;
						}
					}
				}
			}
		}
		for (s_chr = 0; s_chr < ref.size(); s_chr++)
		{
			for (ii = 0; ii < ref[s_chr].loci.size(); ii++)
			{
				for (iii = 0; iii < ref[s_chr].loci[ii].snps.size(); iii++)
				{
					obs_het.chrom_set[s_chr].loci[ii].snp_freq[iii] = obs_het.chrom_set[s_chr].loci[ii].snp_freq[iii] / (num_adults * 2);
				}
			}
		}

		return obs_het;
	}//end calc_obs_het

	Genome Fst(Genome &subset1, Genome &subset2, vector<Chromosome> &ref)
	{

		Genome fst(ref);
		int NumChrom = ref.size();
		double Hs_sub1, Hs_sub2, Ht, Fst, avg_freq;
		int whichchromosome;
		size_t whichallele, marker;

		for (whichchromosome = 0; whichchromosome < NumChrom; whichchromosome++)
		{
			for (marker = 0; marker < ref[whichchromosome].loci.size(); marker++)
			{
				vector <double> averageallelefreqs;
				for (whichallele = 0; whichallele < ref[whichchromosome].loci[marker].snps.size(); whichallele++)
				{
					avg_freq = (subset1.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele] +
						subset2.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele]) / 2;
					averageallelefreqs.push_back(avg_freq);
				}
				Hs_sub1 = 1;
				Hs_sub2 = 1;
				Ht = 1;

				for (whichallele = 0; whichallele < ref[whichchromosome].loci[marker].snps.size(); whichallele++)
				{
					Hs_sub1 = Hs_sub1 - (subset1.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele] *
						subset1.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele]);
					Hs_sub2 = Hs_sub2 - (subset2.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele] *
						subset2.chrom_set[whichchromosome].loci[marker].snp_freq[whichallele]);
					Ht = Ht - (averageallelefreqs[whichallele] * averageallelefreqs[whichallele]);
				}

				if (Ht > 0)
					Fst = 1 - (Hs_sub1 + Hs_sub2) / (2 * Ht);
				else
					Fst = -1.0;

				fst.chrom_set[whichchromosome].loci[marker].snp_freq[0] = Fst;
			}//marker
		}//chromosome
		return fst;
	}//end Fst

	Genome calc_exp_het(vector <Chromosome> &ref, Genome &allele_freq)
	{
		size_t f, ff, al;
		Genome exp_het(ref);
		for (f = 0; f < ref.size(); f++)
		{
			for (ff = 0; ff < ref[f].loci.size(); ff++)
			{
				for (al = 0; al < ref[f].loci[ff].snps.size(); al++)
				{
					exp_het.chrom_set[f].loci[ff].snp_freq.push_back(0);
					//calculate expected heterozygosity
					exp_het.chrom_set[f].loci[ff].snp_freq[al] = exp_het.chrom_set[f].loci[ff].snp_freq[al] +
						(allele_freq.chrom_set[f].loci[ff].snp_freq[al] * allele_freq.chrom_set[f].loci[ff].snp_freq[al]);
				}
				exp_het.chrom_set[f].loci[ff].snp_freq[al] = 1 - exp_het.chrom_set[f].loci[ff].snp_freq[al];
			}
		}
		return exp_het;
	}//end calc heterozygosity

	Genome weighted_fst(Genome &fsts, vector <Chromosome> &ref)
	{
		/***************************************************************************************************************
		Code to weight Fsts to generate smooth genome-wide distributions
		Uses kernel-smoothing average, following the approach in Hohenlohe et al. (2010) and Catchen et al. (2013; Stacks)
		for each region at nucleotide c, contribution of statistic at position p to the region average is
		weighted by Gaussian function:
		exp((-(p-c)*(p-c))/(2*sigma*sigma)), where sigma = 150 kb (chosen arbitrarily)
		the distribution is truncated at 3sigma.
		Shifted the moving average by a step of 100 kb
		Further weight each statistic at each nucletide position by (nk-1) where nk = # alleles sampled at site k
		***************************************************************************************************************/
		size_t chr, loc;
		double p, c, sigma;
		int num_to_avg = 0;
		double dN, avgweights = 0;
		int start, end;
		double coef = 0;
		int locus, y;
		bool calc = false;
		double weightedFst;
		Genome weighted_fsts(ref);

		sigma = 0;
		for (size_t chr = 0; chr < ref.size(); chr++)
		{
			for (loc = 0; loc < ref[chr].loci.size(); loc++)
				sigma++;
		}
		sigma = sigma / 1000;
		locus = 0;
		for (chr = 0; chr < ref.size(); chr++)
		{
			for (loc = 0; loc < ref[chr].loci.size(); loc++)
			{
				if (fsts.chrom_set[chr].loci[loc].snp_freq[0] >= 0)
				{
					//set up weighting parameters
					locus = loc;
					start = locus - 3 * sigma;
					if (start < 0)
						start = 0;
					end = locus + 3 * sigma;
					if (end > ref[chr].loci.size() - 1)
						end = ref[chr].loci.size() - 1;

					//weighted Fst
					num_to_avg = 0;
					avgweights = 0;
					coef = 0;
					for (y = start; y <= end; y++)
					{
						if (fsts.chrom_set[chr].loci[y].snp_freq[0] >= 0)
						{
							p = y;
							c = locus;
							weightedFst = fsts.chrom_set[chr].loci[y].snp_freq[0] * exp(double(-1 * ((p - c)*(p - c))) / (double(2 * sigma*sigma)));
							avgweights = avgweights + weightedFst;
							coef = coef + exp(double((-1 * ((p - c)*(p - c)))) / (double(2 * sigma*sigma)));
							num_to_avg++;
						}
					}
					dN = num_to_avg;
					weightedFst = avgweights / coef;
					weighted_fsts.chrom_set[chr].loci[locus].snp_freq[0] = weightedFst;
				}//end weighting Fsts (if Polymorphisms == 1)
				else
					weighted_fsts.chrom_set[chr].loci[loc].snp_freq[0] = 0;
			}//end x
		}//end w
		return weighted_fsts;
	}//end weighted fst

	Cutoffs calculate_fst_CIs(Genome &fsts, vector <Chromosome> &ref)
	{
		double meanFst, stdevFst, varFst;
		Cutoffs these_cutoffs;
		size_t m, mm;
		int count;
		meanFst = 0;
		count = 0;
		for (m = 0; m < ref.size(); m++)
		{
			for (mm = 0; mm < ref[m].loci.size(); mm++)
			{
				if (fsts.chrom_set[m].loci[mm].snp_freq[0] >= 0)
				{
					meanFst = meanFst + fsts.chrom_set[m].loci[mm].snp_freq[0];
					count++;
				}
			}
		}
		meanFst = meanFst / count;

		varFst = 0;
		count = 0;
		for (m = 0; m < ref.size(); m++)
		{
			for (mm = 0; mm < ref[m].loci.size(); mm++)
			{
				varFst = varFst + (meanFst - fsts.chrom_set[m].loci[mm].snp_freq[0])*(meanFst - fsts.chrom_set[m].loci[mm].snp_freq[0]);
				count++;
			}
		}
		varFst = varFst / count;
		stdevFst = sqrt(varFst);

		these_cutoffs.critical_99 = meanFst + 2.57583 * stdevFst;
		these_cutoffs.critical_98 = meanFst + 2.32635 * stdevFst;
		these_cutoffs.critical_95 = meanFst + 1.95996 * stdevFst;
		these_cutoffs.critical_90 = meanFst + 1.64485 * stdevFst;
		these_cutoffs.critical_80 = meanFst + 1.28155 * stdevFst;

		cout << "99% Critical Value Cutoff is " << these_cutoffs.critical_99 << ".\n";
		return these_cutoffs;
	}//end CalculateFstCIs 

};//end Population