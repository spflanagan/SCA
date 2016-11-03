//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date Updated: 1 Nov 2016
//Date Started: 1 Nov 2016
//Purpose: Simulate mated pairs based on population allele frequencies and compare inferred maternal allele frequencies.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include "random_numbers.h"

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

int main()
{
	int i, ii, iii, count, index, pos, gt_index, line_count, num_offspring, affected_count;
	double error_rate, num_allowed;
	string vcf_name, out_name, line, stemp, dad_allele, mom_allele;
	ifstream vcf;
	ofstream out;
	vector<int> fem_index, mal_index, fem_count, mal_count, inf_count;
	vector<double> fem_af, mal_af, mom_af, inf_af;
	vector<string> refs, alts,ids;
	string chrom, ID, ref, alt, filter, info, format, temp, temp2, last_chr, quality;

	vcf_name = "../../results/stacks/batch_1.vcf";
	out_name = "../../results/sca_simulation_output/maternal_alleles_sim_out_error1.txt";
	num_offspring = 1000;
	error_rate = 1;

	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	line_count = 0;
	while (universal_getline(vcf, line))
	{
		if (!vcf.eof())
		{
			if (line.substr(0, 2) == "#C")
			{
				//need to determine male-offspring indices
				stringstream ss;
				ss.str(line);
				count = index = 0;
				while (ss >> stemp)
				{
					if (count > 8)
					{
						if (stemp.find( "FEM") != stemp.npos)
						{
							fem_index.push_back(index);
						}
						if (stemp.find("PRM") != stemp.npos || stemp.find("NPM") != stemp.npos)
						{
							mal_index.push_back(index);							
						}
						index++;
					}
					count++;
				}//while ss >>stemp
			}
			if (line.substr(0, 1) != "#") //then it's a new locus
			{
				fem_af.push_back(0);
				fem_count.push_back(0);
				mal_af.push_back(0);
				mal_count.push_back(0);
				stringstream ss;
				ss.str(line);
				ss >> chrom >> pos >> ID >> ref >> alt >> quality >> filter >> info >> format;
				refs.push_back(ref);
				alts.push_back(alt);
				ids.push_back(ID);
				stringstream ssi;
				ssi << format;
				count = 0;
				vector<string> format_fields;
				while (getline(ssi, temp, ':'))
				{
					if (temp == "GT") //determine which slot in the FORMAT is the genotype (GT)
						gt_index = count;
					count++;
				}

				index = 0;
				vector<string> alleles;
				while (ss >> stemp)
				{
					stringstream ssii;
					ssii << stemp;
					vector<string> vectemp;
					while (getline(ssii, temp, ':'))
					{
						vectemp.push_back(temp);
					}
					for (i = 0; i < fem_index.size(); i++)
					{
						if (fem_index[i] == index)
						{
							if (vectemp[gt_index] != "./.")
								fem_count[line_count]++;
							if (vectemp[gt_index].substr(0,1) == "0")
							{
								fem_af[line_count]++;
							}
							if (vectemp[gt_index].substr(2,1) == "0")
							{
								fem_af[line_count]++;
							}
						}
					}
					for (i = 0; i < mal_index.size(); i++)
					{
						if (mal_index[i] == index)
						{
							if (vectemp[gt_index] != "./.")
								mal_count[line_count]++;
							if (vectemp[gt_index].substr(0, 1) == "0")
							{
								mal_af[line_count]++;
							}
							if (vectemp[gt_index].substr(2, 1) == "0")
							{
								mal_af[line_count]++;
							}
						}
					}
					index++;
				}
				mal_af[line_count] = mal_af[line_count] / (2*mal_count[line_count]);
				fem_af[line_count] = fem_af[line_count] / (2*fem_count[line_count]);
				line_count++;
			}
		}
	}
	vcf.close();

	cout << "\nFinished parsing vcf file with " << mal_af.size() << " loci.\n";
	num_allowed = error_rate*mal_af.size();
	//generate samples
	for (i = 0; i < mal_af.size(); i++)
	{
		mom_af.push_back(0);
		inf_af.push_back(0);
		inf_count.push_back(0);
	}
	cout << "\nSimulate and infer maternal alleles";
	for (i = 0; i < num_offspring; i++)
	{
		affected_count = 0;
		for (ii = 0; ii < mal_af.size(); ii++)
		{
			//assign parental alleles
			string dad1, dad2,mom1, mom2, kid1, kid2;
			if (genrand() < mal_af[ii])
				dad1 = refs[ii];
			else
				dad1 = alts[ii];
			if (genrand() < mal_af[ii])
				dad2 = refs[ii];
			else
				dad2 = alts[ii];
			if (genrand() < fem_af[ii])
				mom1 = refs[ii];
			else
				mom1 = alts[ii];
			if (genrand() < fem_af[ii])
				mom2 = refs[ii];
			else
				mom2 = alts[ii];
			if (genrand() < 0.5)
				kid1 = dad1;
			else
				kid1 = dad2;
			if (genrand() < 0.5)
				kid2 = mom1;
			else
				kid2 = mom2;
			if (kid2 == refs[ii])
				mom_af[ii]++;
			//incorporate error
			if (kid1 != kid2)
			{
				if (genrand() < error_rate && affected_count < num_allowed)
				{
					if (genrand() < 0.5)
						kid1 = kid2;
					else
						kid2 = kid1;
					affected_count++;
				}
			}
			if (dad1 != dad2)
			{
				if (genrand() < error_rate && affected_count < num_allowed)
				{
					if (genrand() < 0.5)
						dad1 = dad2;
					else
						dad2 = dad1;
					affected_count++;
				}
			}
			//infer maternal alleles
			string mom_allele = ".";
			if (dad1 == dad2 && kid1 == kid2) //the case where both are homozygous
			{
				if (dad1 == kid1)
					mom_allele = dad1;
			}
			if (dad1 == dad2 && kid1 != kid2)//the case where dad is homozygous but kid is het
			{
				if (dad1 == kid1)
					mom_allele = kid2;
				if (dad1 == kid2)
					mom_allele = kid1;
			}
			if (dad1 != dad2 && kid1 == kid2) //the case where dad is het but off is hom
			{
				if (dad1 == kid1 || dad2 == kid1)
					mom_allele = kid1;
			}
			if (dad1 != dad2 && kid1 != kid2)//if they're both hets you can't do anything with it.
				mom_allele = ".";
			if (mom_allele != ".")
			{
				inf_count[ii]++;
				if (mom_allele == refs[ii])
					inf_af[ii]++;
			}
		}//all the loci
		cout << "\nOffspring " << i << " error-ridden alleles: " << affected_count;
	}
	for (i = 0; i < mal_af.size(); i++)
	{
		mom_af[i] = mom_af[i]/num_offspring;
		inf_af[i] = inf_af[i]/inf_count[i];
	}

	//output info
	cout << "\nWriting data to file.\n";
	out.open(out_name);
	out << "Locus\tRef\tAlt\tFemaleAF\tMaleAF\tActualMomAF\tInferredMomAF";
	for (i = 0; i < mal_af.size(); i++)
	{
		out << '\n' <<ids[i] << '\t'<<  refs[i] << '\t' << alts[i] << '\t' << fem_af[i] << '\t' << mal_af[i] << '\t' << mom_af[i] <<'\t' << inf_af[i];
	}
	out.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}