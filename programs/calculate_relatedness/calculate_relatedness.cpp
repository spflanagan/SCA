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
			freq.push_back(0);
		}
		if (found2 == false)
		{
			alleles.push_back(al2);
			freq.push_back(0);
		}
	}

	void caclualte_allele_freqs()
	{
		int j;
		for (j = 0; j < alleles.size(); j++)
			freq[j] = freq[j] / count;
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

	//now compare each individual at each locus
	for (i = 0; i < population.size(); i++)
	{
		for (ii = 0; ii < population.size(); ii++)
		{
			if (i != ii)
			{
				for (iii = 0; iii < reference.size(); iii++)
				{
					
					if (population[i].allele1[iii] == population[ii].allele1[iii])
					{
						//calculate relatedness and likelihood estimators
					}

				}
			}
		}
	}
	cout << "\nDone! Input Integer to Quit.\n";
	cin >> i;
	return 0;
}