//Author: Sarah P. Flanagan
//Purpose: Simulate SNPs and mating to evaluate how many SNPs are necessary and sufficient for parentage

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random_numbers.h"


using namespace std;

class individual
{
public:
	vector<char> maternal;
	vector<char> paternal;
	int mom, dad, num_mates;

	individual()
	{
		maternal = paternal = vector<char>();
		mom = dad = num_mates = int();
	}

	void initialize(int num_snps)
	{
		int j;
		char alleles[4] = { 'A', 'C', 'T', 'G' };

		for (j = 0; j < num_snps; j++)
		{
			maternal.push_back(alleles[randnum(4)]);
			paternal.push_back(alleles[randnum(4)]);
		}
		num_mates = 0;
	}
};

void help_message()
{
	cout << "\n\t\tHelp Menu for ParentageSim\n";
	cout << "\nThis model creates unlinked SNPs for parents and generates a CERVUS input file for their offspring\n";
	cout << "\nFemales mate with one male but males can mate multiply\n";
	cout << "-o:\tOutput file name (simulated_genotypes.txt)\n";
	cout << "-F:\tnumber of females (50)\n";
	cout << "-M:\tnumber of males (50)\n";
	cout << "-S:\tnumber of SNPs (10)\n";
	cout << "-m:\tmaximum number of mates per male (5)\n";
	cout << "-f:\tfecundity, aka number of offspring produced per mating (4)\n";
	cout << "-h:\tPrint this help message.\n";
}

int main(int argc, char*argv[])
{
	int i, ii, iii, num_snps, num_males, num_females, mate, max_num_mates, fecundity, num_encounters, num_offspring;
	vector<individual> males, females;
	string output_name = "simulated_genotypes.txt";
	ofstream output;

	//default settings
	num_snps = 10;
	num_males = 50;
	num_females = 50;
	max_num_mates = 5;
	fecundity = 4;
	//read in parameters
	string tempstring1, tempstring2;
	if (argc == 1)
	{
		help_message();
		return 0;
	}
	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h" || tempstring1 == "--help" || tempstring1 == "--H" || tempstring1 == "--Help" || tempstring1 == "--HELP")
		{
			help_message();
			return 0;
		}
		else
		{			
			for (i = 1; i < argc; i++)
			{
				tempstring1 = argv[i];
				tempstring2 = argv[i+ 1];
				i++;
				if (tempstring1 == "-o")
					output_name = tempstring2;
				if (tempstring1 == "-M")
					num_males = atoi(tempstring2.c_str());
				if (tempstring1 == "-F")
					num_females = atoi(tempstring2.c_str());
				if (tempstring1 == "-m")
					max_num_mates= atoi(tempstring2.c_str());
				if (tempstring1 == "-f")
					fecundity = atoi(tempstring2.c_str());
				if (tempstring1 == "-S")
					num_snps = atoi(tempstring2.c_str());		
			}
		}
	}

	//initialize
	for (i = 0; i < num_males; i++)
	{
		males.push_back(individual());
		males[i].initialize(num_snps);
	}
	for (i = 0; i < num_females; i++)
	{
		females.push_back(individual());
		females[i].initialize(num_snps);
	}

	//mating
	output.open(output_name);
	output << "ID\tMom\tDad";
	for (i = 0; i < num_snps; i++)
		output << "\tA" << i << "\tB" << i;
	num_offspring = 0;
	cout << "\nCommencing mating";
	for (i = 0; i < num_females; i++)
	{
		num_encounters = 0;
		while (num_encounters < 50)
		{
			mate = randnum(num_males);
			if (males[mate].num_mates < max_num_mates)
			{
				for (iii = 0; iii < fecundity; iii++)//females mate with one male
				{									//males can mate multiply
					output << '\n' << "OFFSPRING" << num_offspring << "\tFEM" << i << "\tMAL" << mate;
					for (ii = 0; ii < num_snps; ii++)
					{//should really just output this. 
					 //from the mom
						if (genrand() < 0.5)
							output << '\t' << females[i].maternal[ii];
						else
							output << '\t' << females[i].paternal[ii];
						//from the dad
						if (genrand() < 0.5)
							output << '\t' << males[i].maternal[ii];
						else
							output << '\t' << males[i].paternal[ii];
					}
					num_offspring++;
					if (num_offspring % 100 == 0)
						cout << '\n' << num_offspring << " offspring have been created" << std::flush;
				}
				males[mate].num_mates++;
				females[mate].num_mates++;
			}
			num_encounters++;
		}
	}
	output.close();

	cout << "\nSimulation completed.\n" << std::flush;
	return 0;
}