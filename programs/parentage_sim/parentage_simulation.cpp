//Author: Sarah P. Flanagan
//Purpose: Simulate SNPs and mating to evaluate how many SNPs are necessary and sufficient for parentage

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random_numbers.h"

#include <stdio.h>  /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/ 
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

//from http://www.codebind.com/cpp-tutorial/c-get-current-directory-linuxwindows/
std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

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
	double prop_moms_sampled = 0.05;
	vector<individual> males, females;
	string base_name = "simulated";
	string genotype_name, offspring_name, candmoms_name, canddads_name, crv_name,dir;
	ofstream genotypes,offspring,cand_males,cand_females,crv;

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
					base_name = tempstring2;
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


	dir = GetCurrentWorkingDir();

	//initialize adults
	genotype_name = base_name + "_genotypes.txt";
	genotypes.open(genotype_name);
	genotypes << "ID\tMom\tDad";
	for (i = 0; i < num_snps; i++)
		genotypes << "\tA" << i << "\tB" << i;
	candmoms_name = base_name + "_candidate_mothers.txt";
	canddads_name = base_name + "_candidate_fathers.txt";
	cand_females.open(candmoms_name);
	cand_males.open(canddads_name);
	for (i = 0; i < num_males; i++)
	{
		males.push_back(individual());
		males[i].initialize(num_snps);
		cand_males << "MAL" << i << '\n';
		genotypes << '\n' << "MAL" << i << "\tNA\tNA";
		for (ii = 0; ii < num_snps; ii++)
		{
			genotypes << '\t' << males[i].maternal[ii] << '\t' << males[i].paternal[ii];;
		}
	}
	for (i = 0; i < num_females; i++)
	{
		females.push_back(individual());
		females[i].initialize(num_snps);
		cand_females << "FEM" << i << '\n';
		genotypes << '\n' << "FEM" << i << "\tNA\tNA";
		for (ii = 0; ii < num_snps; ii++) 
		{
			genotypes << '\t' << females[i].maternal[ii] << '\t' << females[i].paternal[ii];;
		}
	}
	

	//mating
	offspring_name = base_name + "_offspring.txt";
	offspring.open(offspring_name);
	offspring << "OffspringID\tMomID\tDadID";	
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
					genotypes << '\n' << "OFFSPRING" << num_offspring << "\tFEM" << i << "\tMAL" << mate;
					offspring << '\n' << "OFFSPRING" << num_offspring << "\tFEM" << i << "\tMAL" << mate;
					
					for (ii = 0; ii < num_snps; ii++)
					{//should really just output this. 
					 //from the mom
						if (genrand() < 0.5)
							genotypes << '\t' << females[i].maternal[ii];
						else
							genotypes << '\t' << females[i].paternal[ii];
						//from the dad
						if (genrand() < 0.5)
							genotypes << '\t' << males[i].maternal[ii];
						else
							genotypes  << '\t' << males[i].paternal[ii];
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
	genotypes.close();
	offspring.close();
	cand_females.close();
	cand_males.close();
	cout << "\nSimulation completed.\n" << std::flush;

	//Output the .crv file
	crv_name = base_name + ".crv";
	crv.open(crv_name);
	crv << "[ProgramInfo]\nProgramName = Cervus\n\nProgramVersion = 3.0\nFileVersion = 3.0.7.0\n";
	crv << "\n[Registration]\nUserName =\nUserCompany =\nCode =\n";
	crv << "\n[FileInfo]\nFileName = " << crv_name << "\nFileType = .crv\nCreationDate = 4 / 26 / 2018 3:48 : 48 PM\n";//I could do this better but it takes work
	crv << "\n[GenotypeFile]\nFileName = " << genotype_name << "\nHeaderRow = 1\nReadLocusNames = 1\nFirstAlleleColumnNumber = 4\nIDColumnNumber = 1\nNLoci = " <<
		num_snps << "\nPropLociTyped = 1\nColumnsPerLocus = 2\nSexColumn = 0\nUnknownSexLabel =\n";
	crv << "\n[CodecFile]\nFileName =\nHeaderRow = 1\nUseSameCodingForAllLoci = 1\nGenotypeFileName = " << genotype_name << '\n';
	crv << "\n[AlleleFrequencySummaryFile]\nFileName = " << dir << base_name << "_afs.txt\nDoHardyWeinberg = 1\nHWMinExpectedFrequency = 5\nUseYatesCorrection = 1\nUseBonferroniCorrection = 1\nDoNullAllele = 1\n";
	crv << "\n[AlleleFrequencyDataFile]\nFileName = " << dir << base_name << ".alf\nHeaderRow = 1\n";
	crv << "\n[SimulationParameters]\nAnalysisType = Maternity\nNOffspring = 10000\nNCandidateMales = 0\nPropCandidateMalesSampled = 0\nNCandidateFemales = " <<
		num_females << "\nPropCandidateFemalesSampled = " << prop_moms_sampled << "\nPropLociTyped = 1\nPropLociMistyped = 0.01\nMinTypedLoci = " <<
		num_snps << "\nCriticalStatisticName = LOD\nTruncateAtZero = 0\nRelaxedConfidence = 90\nStrictConfidence = 95\nSimulateInbreeding = 0\nParentRelatedness = 0" <<
		"\nInbreedingRate = 0\nAlwaysTestSelfing = 0\nSimulateFemaleRelatives = 0\nFemalePropRelatives = 0\nFemaleRelatedTo = Offspring\nFemaleRelatedness = 0\n" <<
		"SimulateMaleRelatives = 0\nMalePropRelatives = 0\nMaleRelatedTo = Offspring\nMaleRelatedness = 0\nUseCorrectedLikelihoods = 1\nUseMistypingRateAsLikelihoodErrorRate = 1\nLikelihoodErrorRate = 0.01\n";
	crv << "\n[SimulationSummaryFile]\nFileName = " << dir << base_name << "_sim.txt\n";
	crv << "\n[SimulationDataFile]\nFileName = " << dir << base_name << ".sim\n";
	crv << "\n[SimulationOutput]\nRepeatSimulation = 0\nNRepeats = 1\nApplyPreviousSimulationData = 0\nGenerateTables = 0\nSaveRawStatisticScores = 0\nGenerateHistograms = 0\n" <<
		"NCategories = 0\nMinStatistic = 0\nMaxStatistic = 0\n";
	crv << "\n[PreviousSimulationDataFile]\nFileName =\n";
	crv << "\n[ParentageParameters]\nAnalysisType = Maternity\nUseSimulationParameters = 1\nCalculateConfidenceLevels = 1\nAlwaysTestSelfing = 0\nMinTypedLoci = " <<
		num_snps << "\nUseCorrectedLikelihoods = 1\nUseMistypingRateAsLikelihoodErrorRate = 0\nLikelihoodErrorRate = 0.01\nCriticalStatisticName = LOD\nTruncateAtZero = 0\n";
	crv << "\n[OffspringFile]\nFileName = " << dir << base_name << "_offspring.txt\nHeaderRow = 1\nOffspringIDColumnNumber = 1\nIncludesKnownParents = 0\nKnownParentIDColumnNumber = 0\n" <<
		"IncludesCandidateParents = 0\nCandidateParentIDColumnNumber = 0\n";
	crv << "\n[CandidateFemaleFile]\nFileName = " << dir << base_name << "_candidate_mothers.txt\nHeaderRow = 0\nCandidateParentFormat = One column for all offspring\n" <<
		"OffspringIDColumnNumber = 0\nCandidateParentIDColumnNumber = 1\n";
	crv << "\n[CandidateMaleFile]\nFileName =\nHeaderRow = 0\nCandidateParentFormat = One column for all offspring\nOffspringIDColumnNumber = 0\nCandidateParentIDColumnNumber = 1\n";
	crv << "\n[ParentageSummaryFile]\nFileName = " << dir << base_name << "_par.txt\n";
	crv << "\n[ParentageDataFile]\nFileName = " << dir << base_name << "_par.csv\nOutputType = The most - likely parent\nSortedBy = Joint LOD score\nIncludeNonExclusionProbabilities = 0\n";
	crv << "\n[IdentityParameters]\nTestSexesSeparately = 0\nMinMatchingLoci = 0\nAllowFuzzyMatching = 0\nMaxMismatchingLoci = 0\nShowAllComparisons = 0\n";
	crv << "\n[IdentitySummaryFile]\nFileName =\n";
	crv << "\n[IdentityDataFile]\nFileName =\n";
	crv << "\n[ConversionSourceFile]\nFileName =\nHeaderRow = 0\nAlleleIDLength = 2\nUseFirstIDAsPopulationName = 0\nPopulationColumnNumber = 0\nIDColumnNumber = 0\n" <<
		"MumIDColumnNumber = 0\nDadIDColumnNumber = 0\nFirstAlleleColumnNumber = 0\nNLoci = 0\nContainsAlleleFrequencyData = 0\n";
	crv << "\n[ConversionSourceAlleleFrequencyDataFile]\nFileName =\nHeaderRow = 0\n";
	crv << "\n[ConversionDestinationFile]\nFileName =\nWriteAlleleFrequencyData = 0\nAlleleIDLength = 2\nWritePopulationNames = 0\n";
	crv.close();
	cout << "\n" << base_name << ".crv file has been created.\n";

	return 0;
}