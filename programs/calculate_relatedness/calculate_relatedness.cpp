//Author: Sarah P. Flanagan
//Date: 7 February 2016
//Purpose: To calculate relatedness among individuals using thousands of haplotypes from next-generation sequencing data

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

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

class relatedness_locus
{
public:
	double rxyl, wl;
};

class locus_info
{
public:
	int count;
	string ID;
	vector<string> alleles;
	vector<double> freq;

	locus_info()
	{
		count = int();
		ID = string();
		alleles = vector<string>();
		freq = vector<double>();
	}

	void update_alleles(string al1, string al2)
	{
		bool found1, found2;
		int i;
		found1 = found2 = false;
		for (i = 0; i < alleles.size(); i++)
		{
			if (alleles[i] == al1)
			{
				freq[i]++;
				found1 = true;
			}
			if (alleles[i] == al2)
			{
				freq[i]++;
				found2 = true;
			}
		}
		if (found1 == false)
		{
			alleles.push_back(al1);
			freq.push_back(1);
		}
		if (found2 == false)
		{
			if (al1 != al2)
			{
				alleles.push_back(al2);
				freq.push_back(1);
			}
			else
				freq.back()++;//you know it was the last one because if al2 wasn't found neither was al1.
		}
	}

	void caclualte_allele_freqs()
	{
		int j;
		for (j = 0; j < alleles.size(); j++)
			freq[j] = freq[j] / (count*2);
	}

	relatedness_locus calc_relatedness(string al_a, string al_b, string al_c, string al_d)
	{
		int iv, Sab, Sbc, Sbd, Sac, Sad;
		double pa, pb;
		relatedness_locus r;
		pa = pb = 0;
		for (iv = 0; iv < alleles.size(); iv++)
		{
			if (alleles[iv] == al_a)
			{
				pa = freq[iv];
			}
			if (alleles[iv] == al_b)
			{
				pb = freq[iv];
			}
		}
		if (al_a == al_b)
			Sab = 1;
		else
			Sab = 0;
		if (al_a == al_c)
			Sac = 1;
		else
			Sac = 0;
		if (al_a == al_d)
			Sad = 1;
		else
			Sad = 0;
		if (al_b == al_c)
			Sbc = 1;
		else
			Sbc = 0;
		if (al_b == al_d)
			Sbd = 1;
		else
			Sbd = 0;
		r.rxyl= ((pa*(Sbc + Sbd)) + (pb*(Sac + Sad)) - (4 * pa*pb)) / ((1 + Sab)*(pa + pb) - (4 * pa*pb));
		r.wl = ((1 + Sab)*(pa + pb) - (4 * pa*pb)) / (2 * pa*pb);
		return r;
	}
};

class individual
{
public:
	string ID;
	vector<string> allele1;
	vector<string> allele2;

	individual()
	{
		ID = string();
		allele1 = vector<string>();
		allele2 = vector<string>();
	}
};

class relatedness_scores
{
public:
	string focal_ID;
	vector<string> comparison_IDs;
	vector<double> r;

	relatedness_scores()
	{
		focal_ID = string();
		comparison_IDs = vector < string > ();
		r = vector<double>();
	}
};
int main()
{
	int i, ii, iii, iv, count, locus_count;
	bool kinship_format;
	string id, allele1, allele2, tmp;
	string kinship_name, allelefreq_name, relatedness_name, line;
	ifstream kinship;
	ofstream allelefreqs, relatedness;
	vector<locus_info> reference;
	vector<individual> population;
	vector<relatedness_scores> r_values;

	relatedness_name = "../../results/relatedness/genotypes99_10loci.rout.txt";
	kinship_name = "../../results/relatedness/genotypes99_10loci.txt";
	kinship_format = false; //if true it's kinship format, if false it's CERVUS format
	
	kinship.open(kinship_name);
	FileTest(kinship, kinship_name);
	count = 0;
	while (universal_getline(kinship, line))
	{
		if (!kinship.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> id;
			if (count == 0)
			{
				if (kinship_format)
				{
					while (ss >> allele1)
					{
						reference.push_back(locus_info());
						reference.back().ID = allele1.substr(0, allele1.size() - 1);
						reference.back().count = 0;
					}
				}
				else
				{
					while (ss >> allele1 >> allele2)
					{
						reference.push_back(locus_info());
						reference.back().ID = allele1.substr(0, allele1.size() - 1);
						reference.back().count = 0;
					}
				}
			}
			else
			{
				population.push_back(individual());
				population.back().ID = id;
				locus_count = 0;
				if (kinship_format)
				{
					while (ss >> allele1)
					{
						stringstream ssa;
						ssa.str(allele1);
						vector<string> all_vec;
						while (getline(ssa, tmp, '/'))
							all_vec.push_back(tmp);
						population.back().allele1.push_back(all_vec[0]);
						population.back().allele2.push_back(all_vec[1]);
						if (all_vec[0] != "0")
						{
							reference[locus_count].update_alleles(all_vec[0], all_vec[1]);
							reference[locus_count].count++;
						}
						locus_count++;
					}
				}
				else
				{
					while (ss >> allele1 >> allele2)
					{
						population.back().allele1.push_back(allele1);
						population.back().allele2.push_back(allele2);
						if (allele1 != "0")
						{
							reference[locus_count].update_alleles(allele1, allele2);
							reference[locus_count].count++;
						}
						locus_count++;
					}
				}//not kinship format
			}
			count++;
		}
	}
	kinship.close();
	cout << "\nSuccessfully read in " << population.size() << " individuals, each with " << reference.size() << " loci.";

	//calculate population allele frequencies
	for (i = 0; i < reference.size(); i++)
		reference[i].caclualte_allele_freqs();

	//set up r_values
	for (i = 0; i < population.size(); i++)
	{
		r_values.push_back(relatedness_scores());
		r_values[i].focal_ID = population[i].ID;
		for (ii = 0; ii < population.size(); ii++)
		{
			if (i != ii)
			{
				r_values[i].r.push_back(0);
				r_values[i].comparison_IDs.push_back(population[ii].ID);
			}
		}
	}

	//now compare each individual at each locus
	relatedness.open(relatedness_name);
	relatedness << "Ind1\tInd2\trxy\tW\tr";
	int index;
	for (i = 0; i < population.size(); i++)
	{
		index = 0;
		for (ii = 0; ii < population.size(); ii++)
		{
			if (i != ii)
			{
				double W, r;
				W = r = 0;
				for (iii = 0; iii < reference.size(); iii++)
				{
					//calculate relatedness where population[i] is the proband
					if (population[i].allele1[iii] != "0" && population[ii].allele1[iii] != "0" && reference[iii].alleles.size() > 1)
					{
						relatedness_locus rxy, ryx;
						rxy = reference[iii].calc_relatedness(population[i].allele1[iii], population[i].allele2[iii], population[ii].allele1[iii], population[ii].allele2[iii]);
						ryx = reference[iii].calc_relatedness(population[ii].allele1[iii], population[ii].allele2[iii], population[i].allele1[iii], population[i].allele2[iii]);
						r = r + (((rxy.rxyl*rxy.wl) + (ryx.rxyl*ryx.wl)) / 2);
						W = W + ((rxy.wl + ryx.wl) / 2);
						//cout << "\n" << reference[iii].ID << '\t' << r;
					}
				}
				if (W > 0)
					r_values[i].r[index] = r / W;
				else
					r_values[i].r[index] = r;
				relatedness << '\n' << population[i].ID << '\t' << population[ii].ID << '\t' << r << '\t' << W << '\t' << r_values[i].r[index];
				index++;
			}
		}
	}
	relatedness.close();
	cout << "\nDone! Input Integer to Quit.\n";
	cin >> i;
	return 0;
}