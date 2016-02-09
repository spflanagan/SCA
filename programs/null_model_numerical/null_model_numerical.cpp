//Author: Sarah P. Flanagan (sflanagan@bio.tamu.edu)
//Date: 13 November 2015
//Purpose: Generate null distributions of Fst based on sampling a single population
//Uses a numerical analysis approach

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <sstream>

using namespace std;
class population
{
public:
	int pop_size, sampled_p;
	double p, q, het;

	population()
	{
		pop_size = int();
		sampled_p = int();
		p = double();
		q = double();
		het = double();
	}
};

class sampled_inds
{
public:
	vector<int> allele1;
	vector<int> allele2;

	sampled_inds()
	{
		allele1 = vector<int>();
		allele2 = vector<int>();
	}
};

class sampled_pop
{
public:
	vector<sampled_inds> inds;

	sampled_pop()
	{
		inds = vector<sampled_inds>();
	}
};

int main(int argc, char* argv[])
{
	int end, time, i;
	int t, tt, pop_size;
	int num_pops, num_reps, num_per_p, num_gens, num_sampled_inds_a, num_sampled_inds_b;
	double migration_rate, pbar, qbar, fst, ht, last_p, Nm;
	double nc, s2, a, b, c, wc_fst, pq, rs, wc_fst8, varp;
	int r, nbar, pbar_index;
	bool interactivemode = false;
	random_device rd;
	default_random_engine generator(rd());
	vector<population> pops;
	population total_pop;
	ofstream output, sample_out;
	string base_name, query, tempstring1, tempstring2;
	stringstream output_name, sample_out_name;

	base_name = "default";
	//num_sampled_inds_a = num_sampled_inds_b = 50;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nnull_model_numerical:\n";
			cout << "Runs numerical analysis to generate a null distribution of Fsts for SCA\n";
			cout << "-o:\tbase file name (include path). Example: a50_b50_\n";
			cout << "-a:\tsample size for subpopulation a\n";
			cout << "-b:\tsample size for subpopulation b\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			cout << "Input integer to quit.\n";
			cin >> end;
			return 0;
		}
		interactivemode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nnull_model_numerical:\n";
			cout << "Runs numerical analysis to generate a null distribution of Fsts for SCA\n";
			cout << "-o:\tbase file name (include path). Example: a50_b50_\n";
			cout << "-a:\tsample size for subpopulation a\n";
			cout << "-b:\tsample size for subpopulation b\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			cout << "Input integer to quit.\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-o")
			base_name = tempstring2;
		if (tempstring1 == "-a")
			num_sampled_inds_a = atoi(tempstring2.c_str());
		if (tempstring1 == "-b")
			num_sampled_inds_b = atoi(tempstring2.c_str());
	}

	if (interactivemode)
	{
		cout << "\nProvide base file name (include path). Example: a50_b50_\n";
		cin >> base_name;
		cout << "\nProvide sample size A\n";
		cin >> num_sampled_inds_a;
		cout << "\nProvide sample size B\n";
		cin >> num_sampled_inds_b;
	}


	num_pops = 1;
	Nm = 1;
	pop_size = 1000;
	total_pop.pop_size = pop_size * num_pops;
	num_reps = 2000;
	num_gens = 5000;
	migration_rate = Nm / pop_size;
	cout << "\nRunning with " << num_pops << " demes, " << " sampling " << num_sampled_inds_a << " and " << num_sampled_inds_b << " from one pop, with migration rate "
		<< migration_rate << " and a population size of " << pop_size << " per deme.\nRunning " << num_reps
		<< " number of reps with " << num_gens << " generations.\n";
	//initialize variables
	for (i = 0; i < num_pops; i++)
	{
		pops.push_back(population());
	}
	sample_out_name << base_name << "sampledpops.txt";
	sample_out.open(sample_out_name.str());
	sample_out << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	output_name << base_name << "output.txt";
	output.open(output_name.str());
	output << "pbar\tHt\tHs\tWrightsFst\tWCFst\tWCFstSimple\tavgp";
	
	population subpop_a, subpop_b;
	double allele1, allele2, total_p, total_het;
	int samp_index = 0;

	for (t = 0; t < num_reps; t++)
	{
		//set migtrant p-value
		uniform_real_distribution <double> unidist(0.05, 0.95);
		pbar = unidist(generator);
		qbar = 1 - pbar;

		for (tt = 0; tt < num_gens; tt++)
		{
			total_pop.p = 0;
			total_pop.q = 0;
			total_pop.het = 0;
			if (tt == 0)
			{
				for (i = 0; i < num_pops; i++)
				{
					pops[i].p = 0.5;
					pops[i].q = 1 - pops[i].p;
					pops[i].het = pops[i].p * pops[i].q * 2;
					pops[i].pop_size = total_pop.pop_size / num_pops;
					total_pop.p = total_pop.p + pops[i].p;
					total_pop.q = total_pop.q + pops[i].q;
					total_pop.het = total_pop.het + pops[i].het;
				}
				total_pop.p = total_pop.p / num_pops;
				total_pop.q = total_pop.q / num_pops;
				total_pop.het = total_pop.het / num_pops;
			}
			nbar = 0;
			for (i = 0; i < num_pops; i++)
			{
				last_p = pops[i].p;
				//drift
				binomial_distribution<int> distribution((2 * pops[i].pop_size), pops[i].p);
				pops[i].sampled_p = distribution(generator);
				pops[i].p = (double)pops[i].sampled_p / (2 * (double)pops[i].pop_size);
				//Island model determines p
				//pt = pt-1(1-m)+pbarm
				pops[i].p = pops[i].p*(1 - migration_rate) + (pbar*migration_rate);
				// don't allow wild changes in allele frequencies of 0.8 or higher
				if (pops[i].p - last_p > 0.8 || pops[i].p - last_p < -0.8)
					pops[i].p = last_p;
				if (pops[i].p < 0)
					pops[i].p = 0;
				if (pops[i].p > 1)
					pops[i].p = 1;
				pops[i].q = 1 - pops[i].p;
				pops[i].het = 2 * pops[i].p * pops[i].q;
				total_pop.p = total_pop.p + pops[i].p;
				total_pop.q = total_pop.q + pops[i].q;
				total_pop.het = total_pop.het + pops[i].het;
				nbar = nbar + pops[i].pop_size;
			}
			total_pop.p = total_pop.p / num_pops;
			total_pop.q = total_pop.q / num_pops;
			total_pop.het = total_pop.het / num_pops;
			nbar = nbar / num_pops;
			//then calculate fst three ways:
			//Wright's: Fst=(Ht-Hs)/Ht aka Nei's: 1-sum(Hs)/HtN
			ht = 1 - ((total_pop.p*total_pop.p) + (total_pop.q*total_pop.q));
			if (ht > 0)
				fst = (ht - total_pop.het) / ht;
			else
				fst = 0;
			//fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
			//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
			//Fst=a/(a+b+c)
			//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
			//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
			//c=0.25hbar
			//nc=(rnbar-sum(n2/rnbar))/(r-1)
			r = num_pops;
			nc = r*nbar;
			s2 = varp = 0;
			for (i = 0; i < num_pops; i++)
			{
				s2 = s2 + pops[i].pop_size*(pops[i].p - total_pop.p)*(pops[i].p - total_pop.p) / ((r - 1)*nbar);
				nc = nc - ((pops[i].pop_size*pops[i].pop_size) / (r*nbar));
				varp = varp + (total_pop.p - pops[i].p)*(total_pop.p - pops[i].p);
			}
			varp = varp / num_pops;
			wc_fst8 = varp / (pbar*(1 - pbar));
			nc = nc / (r - 1);
			pq = (total_pop.p*(1 - total_pop.p));
			rs = ((r - 1) / r)*s2;
			a = (nbar / nc)*(s2 - (1 / (nbar - 1))*(pq - rs - (0.25*total_pop.het)));
			b = (nbar / (nbar - 1))*(pq - rs - ((2 * nbar - 1) / (4 * nbar))*total_pop.het);
			c = 0.5*total_pop.het;
			if ((a + b + c)>0)
				wc_fst = a / (a + b + c);
			else
				wc_fst = 0;
		}//gens
		output << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;

		//sample populations
		nbar = 0;
		//choose a pop
		uniform_real_distribution <double> popdist(0, num_pops);
		samp_index = int(popdist(generator));
		total_pop.het = total_pop.p = total_pop.q = 0;
		subpop_a.pop_size = pops[samp_index].pop_size;
		nbar = nbar + subpop_a.pop_size;
		subpop_a.p = 0;	
		for (tt = 0; tt < num_sampled_inds_a; tt++)
		{
			uniform_real_distribution <double> alleledist(0, 1);
			allele1 = alleledist(generator);
			allele2 = alleledist(generator);
			if (allele1 <= pops[samp_index].p)
				subpop_a.p++;
			if (allele2 <= pops[samp_index].p)
				subpop_a.p++;
		}
		subpop_a.p = subpop_a.p / (2 * num_sampled_inds_a);
		subpop_a.q = 1 - subpop_a.p;
		subpop_a.het = 2 * subpop_a.p * subpop_a.q;
		total_pop.p = total_pop.p + subpop_a.p;
		total_pop.q = total_pop.q + subpop_a.q;
		total_pop.het = total_pop.het + subpop_a.het;
	
		subpop_b.pop_size = pops[samp_index].pop_size;
		nbar = nbar + subpop_b.pop_size;
		subpop_b.p = 0;
		for (tt = 0; tt < num_sampled_inds_b; tt++)
		{
			uniform_real_distribution <double> unidist(0, 1);
			allele1 = unidist(generator);
			allele2 = unidist(generator);
			if (allele1 <= pops[samp_index].p)
				subpop_b.p++;
			if (allele2 <= pops[samp_index].p)
				subpop_b.p++;
		}
		subpop_b.p = subpop_b.p / (2 * num_sampled_inds_b);
		subpop_b.q = 1 - subpop_b.p;
		subpop_b.het = 2 * subpop_b.p * subpop_b.q;
		total_pop.p = total_pop.p + subpop_b.p;
		total_pop.q = total_pop.q + subpop_b.q;
		total_pop.het = total_pop.het + subpop_b.het;
			
			
		total_pop.p = total_pop.p / 2;
		total_pop.q = total_pop.q / 2;
		total_pop.het = total_pop.het / 2;
		ht = 1 - ((total_pop.p*total_pop.p) + (total_pop.q*total_pop.q));
		if (ht > 0)
			fst = (ht - total_pop.het) / ht;
		else
			fst = 0;
		//fst = (roundf(fst * 1000000) / 1000000);//keep 6 decimal points
		//Weir and Cockerham: Fst = (f0-f1)/(1-f1) where 1-f1 is heterozygosity
		//Fst=a/(a+b+c)
		//a=(nbar/nc)(s2-(1/(nbar-1))(pbar(1-pbar)-((r-1)/r)*s2-0.25hbar))
		//b=(nbar/(nbar-1))*(pbar*(1-pbar)-((r-1)/r)*s2-((2nbar-1)/4nbar)*hbar)
		//c=0.25hbar
		//nc=(rnbar-sum(n2/rnbar))/(r-1)
		r = 2;
		nbar = nbar / 2;
		nc = r*nbar;
		s2 = varp = 0;
		s2 = s2 + subpop_a.pop_size*(subpop_a.p - total_pop.p)*(subpop_a.p - total_pop.p) / ((r - 1)*nbar);
		nc = nc - ((subpop_a.pop_size*subpop_a.pop_size) / (r*nbar));
		varp = varp + (total_pop.p - subpop_a.p)*(total_pop.p - subpop_a.p);
		s2 = s2 + subpop_b.pop_size*(subpop_b.p - total_pop.p)*(subpop_b.p - total_pop.p) / ((r - 1)*nbar);
		nc = nc - ((subpop_b.pop_size*subpop_b.pop_size) / (r*nbar));
		varp = varp + (total_pop.p - subpop_b.p)*(total_pop.p - subpop_b.p);

		varp = varp / r;
		wc_fst8 = varp / (pbar*(1 - pbar));
		nc = nc / (r - 1);
		pq = (total_pop.p*(1 - total_pop.p));
		rs = ((r - 1) / r)*s2;
		a = (nbar / nc)*(s2 - (1 / (nbar - 1))*(pq - rs - (0.25*total_pop.het)));
		b = (nbar / (nbar - 1))*(pq - rs - ((2 * nbar - 1) / (4 * nbar))*total_pop.het);
		c = 0.5*total_pop.het;
		if ((a + b + c)>0)
			wc_fst = a / (a + b + c);
		else
			wc_fst = 0;
		sample_out << '\n' << pbar << '\t' << ht << '\t' << total_pop.het << '\t' << fst << '\t' << wc_fst << '\t' << wc_fst8 << '\t' << total_pop.p;
	}//reps
	output.close();
	sample_out.close();



	if (interactivemode)
	{
		cout << "\nDone! Input integer to quit.\n";
		cin >> end;
		return 0;
	}
	else
	{
		cout << "\n";
		return 0;
	}
}