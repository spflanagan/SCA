//Author: Sarah P. Flanagan
//Date: 15 February 2016
//Purpose: The purpose of this program is to calculate relatedness as "band-sharing",
//aka the number of shared alleles between a known parent (father) and offspring, divided by 2N (N = number of loci)
//While calculating band-sharing, this program also outputs the number of incompatible loci between fathers and offspring
//Input: one file containing father-offspring pairs, one cervus-format file

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

string find_and_replace(string &s, string toReplace, string replaceWith)
{
	size_t pos = 0;
	while ((pos = s.find(toReplace, pos)) != string::npos) {
		s.replace(pos, toReplace.length(), replaceWith);
		pos += replaceWith.length();
	}
	return s;
}

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

class par_off_pair
{
public:
	int off_index, par_index, num_loci;
	double shared, incompatible;
	string parent, offspring;

	par_off_pair()
	{
		off_index = par_index = num_loci = int();
		shared = incompatible = double();
		parent = offspring = string();
	}
};

int calc_shared(string al_a, string al_b, string al_c, string al_d)
{
	int iv, Sab, Sbc, Sbd, Sac, Sad;
	int shared;
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
	shared = (Sac + Sad + Sbc + Sbd);
	return shared;
}

int main()
{
	int i, ii,iii, index, count, locus_count;
	string allele1, allele2, id, chrom, pos;
	string kinship_name, output_name, pairs_name, line, stemp;
	ifstream kinship, pairs_file;
	ofstream output;
	vector<individual> population;
	vector<par_off_pair> pairs;
	bool found,vcf;
	string empty_char;
	
	pairs_name = "../../results/relatedness/pairwise.combinations.txt";
	output_name = "../../results/relatedness/dradPrunedHaps.allcombos.bandsharing.txt";
	kinship_name = "../../results/parentage/dradPrunedHaps.txt";
	vcf = false;

	if (vcf)
		empty_char = ".";
	else
		empty_char = "0";

	pairs_file.open(pairs_name);
	FileTest(pairs_file, pairs_name);
	while (universal_getline(pairs_file, line))
	{
		if (!pairs_file.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> allele1 >> allele2;
			pairs.push_back(par_off_pair());
			pairs.back().parent = allele1;
			pairs.back().offspring = allele2; 
			pairs.back().incompatible = 0;
			pairs.back().shared = 0;
		}
	}
	pairs_file.close();

	
	kinship.open(kinship_name);
	FileTest(kinship, kinship_name);
	count = locus_count = 0;
	while (universal_getline(kinship, line))
	{
		if (!kinship.eof())
		{
			if (vcf)
			{
				if (line.substr(0, 5) == "CHROM")
				{
					//need to determine male-offspring indices
					stringstream ss;
					ss.str(line);
					count = index = 0;
					while (ss >> stemp)
					{
						if (count > 1)
						{//need to identify which individual is where

							population.push_back(individual());
							population.back().ID = stemp;

							for (i = 0; i < pairs.size(); i++)
							{
								if (stemp == pairs[i].parent)
									pairs[i].par_index = count;//then I can just use ind[indices[count]] to access the individuals if necessary.
									
								if (stemp == pairs[i].offspring)
									pairs[i].off_index = count;//then I can just use ind[indices[count]] to access the individuals if necessary.
									
							}
						}
						count++;
					}//while ss >>stemp
					locus_count = 0;
				}
				else//then it's a new locus
				{

					stringstream ss;
					ss.str(line);
					ss >> chrom >> pos;
					index = 0;
					vector<string> alleles;
					while (ss >> stemp)
					{
						alleles.push_back(stemp);
						find_and_replace(stemp, "/", "\n");
						stringstream ssa;
						ssa.str(stemp);
						getline(ssa, stemp);
						allele1 = stemp;
						getline(ssa, stemp);
						allele2 = stemp;

						population[index].allele1.push_back(allele1);
						population[index].allele2.push_back(allele2);
						index++;
					}//end of reading in individuals
					
					locus_count++;
				}
			}//if vcf
			else
			{
				stringstream ss;
				ss.str(line);
				ss >> id;
				if (count != 0)
				{
					found = false;
					for (i = 0; i < pairs.size(); i++)
					{
						if (id == pairs[i].parent)
						{
							found = true;
							pairs[i].par_index = population.size();
						}
						if (id == pairs[i].offspring)
						{
							found = true;
							pairs[i].off_index = population.size();
						}
					}
					if (found)
					{//only store the relevant individuals
						population.push_back(individual());
						population.back().ID = id;
						locus_count = 0;

						while (ss >> allele1 >> allele2)
						{
							population.back().allele1.push_back(allele1);
							population.back().allele2.push_back(allele2);
							locus_count++;
						}
					}
				}
				count++;
			}//not vcf
			
			if (locus_count == 1000)
			{
				cout << locus_count << " loci successfully parsed.\n";
				locus_count = population[0].allele1.size();
				for (i = 0; i < pairs.size(); i++)
				{
					for (ii = 0; ii < locus_count; ii++)
					{
						if (population[pairs[i].par_index].allele1[ii] != empty_char)
						{
							int sum = 0;
							sum = calc_shared(population[pairs[i].off_index].allele1[ii], population[pairs[i].off_index].allele1[ii],
								population[pairs[i].par_index].allele1[ii], population[pairs[i].par_index].allele2[ii]);
							pairs[i].shared = pairs[i].shared + sum;
							pairs[i].num_loci++;
							//incompatible loci
							if (sum == 0)
								pairs[i].incompatible++;

						}
					}
				}//pairs
				for (i = 0; i < population.size(); i++)
				{
					population[i].allele1.resize(0);
					population[i].allele2.resize(0);
				}
				locus_count= 0;
			}//locus_count
		}//kinship
	}
	kinship.close();
	cout << "\nWriting Data to file.";
	output.open(output_name);
	output << "Father\tOffspring\tNumLoci\tShared\tIncompatible";
	locus_count = population[0].allele1.size();
	for (i = 0; i < pairs.size(); i++)
	{
		for (ii = 0; ii < locus_count; ii++)
		{
			if (population[pairs[i].par_index].allele1[ii] != empty_char)
			{
				int one, two;//offspring alleles one and two
				one = two = 0;//we want to make sure we don't count the same allele twice.
				//shared loci
				if (population[pairs[i].off_index].allele1[ii] == population[pairs[i].par_index].allele1[ii])
				{
					one = 1;
					pairs[i].shared++;
				}
				if (one == 0 && population[pairs[i].off_index].allele1[ii] == population[pairs[i].par_index].allele2[ii])
				{
					pairs[i].shared++;
					one = 2;
				}
				if (one != 1 && population[pairs[i].off_index].allele2[ii] == population[pairs[i].par_index].allele1[ii])
				{
					pairs[i].shared++;
					two = 1;
				}
				if (one != 2 && two == 0 && population[pairs[i].off_index].allele2[ii] == population[pairs[i].par_index].allele2[ii])
				{
					pairs[i].shared++;
					two = 2;
				}

				pairs[i].num_loci++;
				//incompatible loci
				if (one == 0 && two == 0)
					pairs[i].incompatible++;

			}
		}
		pairs[i].shared = pairs[i].shared / (pairs[i].num_loci * 2);//2 times the number of bands that two fingerprints have in common, divided by the total number of bands that the two genotypes have (or #bands/2#loci)
		pairs[i].incompatible = pairs[i].incompatible / pairs[i].num_loci;//the proportion of loci that are incompatible.
		output << '\n' << pairs[i].parent << '\t' << pairs[i].offspring << '\t' << pairs[i].num_loci << '\t' << pairs[i].shared << '\t' << pairs[i].incompatible;
	}//pairs
	output.close();
	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}