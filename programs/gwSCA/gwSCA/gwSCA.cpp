//Date: 14 May 2015
//Adapted from sparrow_project
//Author: Sarah P. Flanagan, sflanagan@bio.tamu.edu
//Purpose: to use a .ped, .map, and an individual data file to:
//(1) sort individuals into groups
//(2) calculate pairwise Fsts using SNP data between the groups
//(3) identify outlier loci using a genome-wide approach.
//Usage:
//sparrow_analysis
//analyzes data in PLINK format.
//-m: map file
//-p: ped file
//-i: individuals data file
//-o: output file name (include path)
//-c: cutoff threshold desired [99, 98, 95, 90, or 80]
//-h: display this message
//no arguments: interactive mode

//this program can either be run in interactive mode (responding to prompts) or by running it with the appropriate flags.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "gwSCA_classes.h"

using namespace std;

int main(int argc, char* argv[])
{
	int end, count, j, snp_count, num_pops, tot_num_inds;
	size_t i, ii, iii, iv, v, s_chr;
	vector <string> map, ped;
	string  line, base_out_name, path_name, out_name;
	stringstream map_name, ped_name, ind_name, ind_info_out, reference_out, summary_name, allele_name;
	ifstream map_file, ped_file, ind_file;
	ofstream ind_info_out_file, ref_file, summary_output, allele_out;
	int yr, loc, sex, age, phen1, phen2, phen3, pop_id, maj_all_loc_1, maj_all_loc_2;
	double gen_dist, maj_all_freq_1, maj_all_freq_2;
	int bp, curr_loc;
	string chr, ID, pop, line2, snp_id, ind_id, last_chr;
	istringstream lin, lin2;
	vector <Genome> pop_afreqs_all, obs_hets;
	vector <Population> pops;
	vector <string> pop_names;
	vector <Chromosome> reference;
	int phen_of_interest;
	vector <Genome> exp_het_all;
	vector < vector <Genome> > pop_afreqs_sel, obs_hets_s, exp_het_sel;
	vector <int> ind_count;
	bool pop_id_found, found;
	int fam_id, pat_id, mat_id, phen, curr_ind;
	char allele_in;
	int chr_index, stringency_value;
	int threshold;//phenotype threshold

	//path_name = "default";// "C:\\Users\\sflanagan\\Dropbox\\sparrows\\";
	//map_name << path_name << "default"; //"All_Inds_Genotyped_On_10K_Chip.map"; //"test_data.map";//
	//ped_name << path_name << "default"; // "All_Inds_Genotyped_On_10K_Chip.ped"; //"test_data.ped";//
	ind_name << "default"; // "SNP-Genotypes-n2348-AlHeLeVeLDstudy-FileToSarahFlanagan.csv"; //"test_data.ind.txt";//
	base_out_name = "gw_sca_";
	phen_of_interest = 4;//add to input
	stringency_value = 99;
	threshold = 0;
	bool phen1_type = true;

	bool interactivemode = true;

	string query;
	string tempstring1, tempstring2;
	map_name << "default";
	ped_name << "default";

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nsparrow analysis:\n";
			cout << "calculates fsts from PLINK data\n";
			cout << "-m:\tmap file (include path)\n";
			cout << "-p:\tped file (include path)\n";
			cout << "-i:\tindividual info file (include path)\n";
			cout << "-o:\toutput file name (include path)\n";
			cout << "-c:\tcutoff threshold desired [99, 98, 95, 90, or 80]\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nsparrow analysis:\n";
			cout << "calculates fsts from PLINK data\n";
			cout << "-m:\tmap file (include path)\n";
			cout << "-p:\tped file (include path)\n";
			cout << "-i:\tindividual info file (include path)\n";
			cout << "-o:\toutput file name (include path)\n";
			cout << "-c:\tcutoff threshold desired [99, 98, 95, 90, or 80]\n";
			cout << "-h:\tdisplay this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-m")
			map_name.str(tempstring2);
		if (tempstring1 == "-p")
			ped_name.str(tempstring2);
		if (tempstring1 == "-i")
			ind_name.str(tempstring2);
		if (tempstring1 == "-o")
			base_out_name = tempstring2;
		if (tempstring1 == "-c")
			stringency_value = atoi(tempstring2.c_str());
	}

	if (interactivemode)
	{
		cout << "\nInput file name of the map file (including path):\n";
		cin >> tempstring1;
		map_name.str(tempstring1);
		cout << "\nInput file name of the ped file (including path):\n";
		cin >> tempstring1;
		ped_name.str(tempstring1);
		cout << "\nInput file name of the individual info file (including path):\n";
		cin >> tempstring1;
		ind_name.str(tempstring1);
		cout << "\nDefault output base name is gw_sca_*. Do you want to change it? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" || tempstring1 == "y")
		{
			cout << "\nInput alternative output base name (including path):\n";
			cin >> base_out_name;
		}
		cout << "\nDefault Fst outlier cutoff threshold is 99%. Do you want to change it? Y or N\n";
		cin >> tempstring1;
		if (tempstring1 == "Y" || tempstring1 == "y")
		{
			cout << "\nInput preferred Fst outlier cutoff threshold:\n";
			cin >> tempstring2;
			stringency_value = atoi(tempstring2.c_str());
		}

		cout << "\n\nmap:\t" << map_name.str();
		cout << "\nped:\t" << ped_name.str();
		cout << "\nindividual info:\t" << ind_name.str();
		cout << "\noutput name:\t" << base_out_name;

		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;
		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}
	
	if (map_name.str() == "default")
	{
		cout << "file name of the map file (including path):\n";
		cin >> tempstring1;
		map_name.str(tempstring1);
		interactivemode = true;
	}

	if (ped_name.str() == "default")
	{
		ped_name.clear();
		cout << "\nfile name of the ped file (including path):\n";
		cin >> tempstring1;
		ped_name.str(tempstring1);
		interactivemode = true;
	}
	if (ind_name.str() == "default")
	{
		ind_name.clear();
		cout << "\nfile name of the individual info file (including path):\n";
		cin >> tempstring1;
		ind_name.str(tempstring1);
		interactivemode = true;
	}

	cout << "\n\nmap:\t" << map_name.str();
	cout << "\nped:\t" << ped_name.str();
	cout << "\nindividual info:\t" << ind_name.str();
	cout << "\noutput name:\t" << base_out_name;
	cout << "\n\nProceeding...\n";

	ind_info_out << base_out_name << "ind_info.txt";
	reference_out << base_out_name << "reference_info.txt";
	summary_name << base_out_name << "summary_output.txt";
	allele_name << base_out_name << "alleles.txt";



	//this is where the stuff happens
	ind_file.open(ind_name.str());
	FileTest(ind_file, ind_name.str());
	count = 0;
	pop_id = 0;
	num_pops = 0;
	tot_num_inds = 0;
	cout << "Parsing Individual Information File.\n";
	while (universal_getline(ind_file, line))
	{
		line2 = find_and_replace(line, ",", "\t");
		lin.clear();
		lin.str(line2);
		if (count == 0)
			lin >> line2;
		else
		{
			if (lin >> ID >> yr >> loc >> pop >> sex >> age >> phen1 >> phen2 >> phen3)
			{
				if (tot_num_inds == 0)
				{
					pops.push_back(Population());
					pop_id = 0;
					pops[0].name = pop;
					pop_names.push_back(pop);
					ind_count.push_back(0);
					num_pops++;
					curr_ind = 0;
				}
				pop_id_found = false;
				for (i = 0; i < pop_names.size(); i++)
				{
					if (pop == pop_names[i])
					{
						pop_id = i;
						pop_id_found = true;
					}
					if (i == pop_names.size() - 1 && pop_id_found == false)
					{
						pop_names.push_back(pop);
						pops.push_back(Population());
						pops[pop_names.size() - 1].name = pop;
						pop_id = pop_names.size() - 1;
						pop_id_found = true;
						ind_count.push_back(0);
						num_pops++;
					}
				}
				//cout << pop_id << '\t';
				//need to figure out if this adult already exists in the database
				if (ind_count[pop_id] == 0)
				{
					curr_ind = 0;
					pops[pop_id].adults[curr_ind].ID = ID;
					pops[pop_id].adults[curr_ind].year[0] = yr;
					if (sex == 1)
						pops[pop_id].adults[curr_ind].female = false;
					else
						pops[pop_id].adults[curr_ind].female = true;
					ind_count[pop_id]++;
					tot_num_inds++;
				}
				else
				{
					curr_ind = -5;
					for (i = 0; i < ind_count[pop_id]; i++)
					{
						if (pops[pop_id].adults[i].ID == ID)
						{
							curr_ind = i;
						}
					}
					if (curr_ind == -5)
					{
						curr_ind = ind_count[pop_id];
						pops[pop_id].adults.push_back(Individual());
						pops[pop_id].adults[curr_ind].ID = ID;
						pops[pop_id].adults[curr_ind].year[0] = yr;
						pops[pop_id].adults[curr_ind].phenotype2.push_back(phen2);
						if (phen1 == 1)//this one remained on the site
							pops[pop_id].adults[curr_ind].phenotype1[0] = 1;
						else
							pops[pop_id].adults[curr_ind].phenotype1[0] = 0;
						if (phen3 == 1)
							pops[pop_id].adults[curr_ind].phenotype3[0] = 1;
						else
							pops[pop_id].adults[curr_ind].phenotype3[0] = 0;
						if (sex == 1)
							pops[pop_id].adults[curr_ind].female = false;
						else
							pops[pop_id].adults[curr_ind].female = true;
						ind_count[pop_id]++;
						tot_num_inds++;
					}//end if(curr_ind == -5
					else
					{
						pops[pop_id].adults[curr_ind].year.push_back(yr);
						pops[pop_id].adults[curr_ind].phenotype2.push_back(phen2);
						if (phen1 == 1)//this one remained on the site
							pops[pop_id].adults[curr_ind].phenotype1.push_back(1);
						else
							pops[pop_id].adults[curr_ind].phenotype1.push_back(0);
						if (phen3 == 1)
							pops[pop_id].adults[curr_ind].phenotype3.push_back(1);
						else
							pops[pop_id].adults[curr_ind].phenotype3.push_back(0);
					}
				}

				
			}//end if...read in line
			else
			{
				cout << "WARNING: Line " << count << " unable to be read:\n" << line << '\n' << line2 << '\n';
			}
		}
		count++;
	}//end while
	ind_file.close();
	cout << "Individual information successfully read.\n";
	cout << "Read data from " << num_pops << " populations: ";
	for (j = 0; j < pop_names.size(); j++)
		cout << pop_names[j] << ", ";
	cout << " and a total of " << tot_num_inds << " individuals.\n\n";

	map_file.open(map_name.str());
	FileTest(map_file, map_name.str());
	count = 0;
	last_chr = -5;
	chr_index = 0;
	snp_count = 0;
	cout << "Parsing .map file.\n";
	while (universal_getline(map_file, line))
	{
		lin.clear();
		lin.str(line);

		if (lin >> chr >> snp_id >> gen_dist >> bp)
		{
			if (snp_count == 0)
			{
				//cout << count << '\t';
				reference.push_back(Chromosome());
				reference[chr_index].ID = chr;
				reference[chr_index].loci[count].snp_name = snp_id;
				reference[chr_index].loci[count].distance = gen_dist;
				last_chr = chr;
				snp_count++;
				count++;
				//	cout << chr << '\t' << chr_index << '\t' << last_chr << '\t' << count << '\n';
			}
			else
			{
				if (chr != last_chr)
				{
					//cout << count << '\t';
					reference.push_back(Chromosome());
					chr_index++;
					count = 0;
					last_chr = chr;
					reference[chr_index].ID = chr;
					reference[chr_index].loci[count].snp_name = snp_id;
					reference[chr_index].loci[count].distance = gen_dist;
					snp_count++;
					count++;
					//if (chr_index > last_chr + 1 && last_chr > 0)
					//{
					//	cout << chr << '\t' << chr_index << '\t' << last_chr << '\t' << count << '\n';
					//	//cout << "\nskipped chromosome " << last_chr + 1 + initial_chr << " and went to " << chr_index + initial_chr << '\n';
					//	reference.push_back(Chromosome());
					//	count = 0;
					//	cout << "Difference: " << chr_index - last_chr << '\n';
					//	initial_chr = initial_chr + (chr_index - last_chr) - 1;
					//	cout << "initial_chr changed to " << initial_chr << '\n';
					//	chr_index = chr - initial_chr;
					//}
					//else
					//{
					//	cout << chr << '\t' << chr_index << '\t' << last_chr << '\t' << count << '\n';
					//	reference.push_back(Chromosome());
					//	count = 0;
					//	chr_index = chr - initial_chr;
					//	last_chr = chr_index;
					//	reference[chr_index].ID = chr;
					//	reference[chr_index].loci[count].snp_name = snp_id;
					//}
					//count++;
					//snp_count++;

				}
				else
				{
					//chr_index = chr - initial_chr;
					reference[chr_index].loci.push_back(Locus_info());
					reference[chr_index].loci[count].snp_name = snp_id;
					reference[chr_index].loci[count].distance = gen_dist;
					count++;
					snp_count++;
				}
			}
		}
	}//end getline(map_file)
	//	cout << count << '\n';
	map_file.close();
	cout << ".map File Successfully Read.\n";

	//cout << "reference size: " << reference.size() << '\n';


	cout << "Read in data for " << snp_count << " SNPs from " << reference.size() << " chromosomes.\n\n";
	//cout << "reference size: " << reference.size() << '\n';
	for (i = 0; i < pops.size(); i++)
	{
		for (ii = 0; ii < pops[i].adults.size(); ii++)
		{
			for (iii = 0; iii < reference.size(); iii++)
			{
				if (reference[iii].loci.size() > 0)
				{
					pops[i].adults[ii].allele1.push_back(Chromosome());
					pops[i].adults[ii].allele2.push_back(Chromosome());
					for (iv = 1; iv < reference[iii].loci.size(); iv++)//first one already initialized when Chromosome() pushed back
					{
						reference[iii].loci[iv].snps.resize(1, 'n');
						reference[iii].loci[iv].snp_freq.resize(1, 0);
						pops[i].adults[ii].allele1[iii].loci.push_back(Locus_info());
						pops[i].adults[ii].allele2[iii].loci.push_back(Locus_info());
					}
				}
			}
		}
	}

	ped_file.open(ped_name.str());
	FileTest(ped_file, ped_name.str());
	cout << "Reading .ped File\n";
	count = 0;
	while (universal_getline(ped_file, line))
	{
		line2 = find_and_replace(line, " ", "\t");
		lin.clear();
		lin.str(line2);
		lin >> fam_id >> ind_id >> pat_id >> mat_id >> sex >> phen;
		found = false;
		while (found == false)
		{
			for (i = 0; i < pops.size(); i++)
			{
				for (ii = 0; ii < pops[i].adults.size(); ii++)
				{
					if (pops[i].adults[ii].ID == ind_id)
					{
						found = true;
						curr_ind = ii;
						pop_id = i;
						//cout << "matched individual " << pops[i].adults[ii].ID << " in pop " << pop_id << '\n';
					}
					if (found == false)
					{
						found = true;
						curr_ind = -5;
					}
				}//end of finding this pop and individual
			}
		}
		count = 0;
		if (curr_ind >= 0)
		{
			for (i = 0; i < reference.size(); i++)
			{
				for (ii = 0; ii < reference[i].loci.size(); ii++)
				{
					lin >> allele_in;
					pops[pop_id].adults[curr_ind].allele1[i].loci[ii].snps[0] = allele_in;
					reference[i].update_chromosome(allele_in, ii);

					lin >> allele_in;
					pops[pop_id].adults[curr_ind].allele2[i].loci[ii].snps[0] = allele_in;
					reference[i].update_chromosome(allele_in, ii);
				}
			}
		}
	}
	ped_file.close();
	//for (i = 0; i < pops.size(); i++)
	//{
	//	pops[i].remove_defaults();
	//	for (ii = 0; ii < pops[i].adults.size(); ii++)
	//	{
	//		for (iii = 0; iii < reference.size(); iii++)
	//		{
	//			for (iv = 0; iv < reference[iii].loci.size(); iv++)
	//			{
	//				cout << pops[i].adults[ii].allele1[iii].loci[iv].snps[0] << '\t';
	//				cout << pops[i].adults[ii].allele2[iii].loci[iv].snps[0] << '\t';
	//			}
	//		}
	//	}
	//}
	for (i = 0; i < reference.size(); i++)
		reference[i].remove_defaults();
	cout << ".ped file has been processed. Writing individual data to file.\n";

	//check if the data are all there
	ind_info_out_file.open(ind_info_out.str());
	if (ind_info_out_file.is_open())
		cout << ind_info_out.str() << " is open!\n";
	ind_info_out_file << "Pop\tIndID\tSex\tNumYears";
	for (i = 0; i < reference.size(); i++)
	{
		for (ii = 0; ii < reference[i].loci.size(); ii++)
		{
			ind_info_out_file << '\t' << reference[i].ID << "." << reference[i].loci[ii].snp_name;
		}
	}
	ind_info_out_file << '\n';
	for (i = 0; i < pops.size(); i++)
	{
		for (ii = 0; ii < pops[i].adults.size(); ii++)
		{
			ind_info_out_file << pops[i].name << '\t' << pops[i].adults[ii].ID << '\t';
			if (pops[i].adults[ii].female == true)
				ind_info_out_file << "F\t";
			else
				ind_info_out_file << "M\t";
			ind_info_out_file << pops[i].adults[ii].year.size();
			for (s_chr = 0; s_chr < reference.size(); s_chr++)
			{
				//  cout << pops[i].adults[ii].allele1[chr].snps.size() << '\t';
				for (iii = 0; iii < pops[i].adults[ii].allele1[s_chr].loci.size(); iii++)
				{
					if (pops[i].adults[ii].allele1[s_chr].loci[iii].snps.size() > 0){
						ind_info_out_file << '\t' << pops[i].adults[ii].allele1[s_chr].loci[iii].snps[0] << "/" << pops[i].adults[ii].allele2[s_chr].loci[iii].snps[0];
						//cout << pops[i].adults[ii].allele1[s_chr].loci[iii].snps[0] << '\t';
					}
					/*if(pops[i].adults[ii].allele1[s_chr].loci[iii].snps[0] != 'n')
					cout << "pop " << pops[i].name << ", adult " << pops[i].adults[ii].ID << ", chrom " << s_chr << "locus " << iii << "was matched!\n";*/
				}
			}
			ind_info_out_file << '\n';
		}
	}
	ind_info_out_file.close();
	cout << "Individual info file written.\n";

	//for (i = 0; i < reference.size(); i++)
	//{
	//	for (ii = 0; ii < reference[i].loci.size(); ii++)
	//	{
	//		for (iii = 0; iii < reference[i].loci[ii].snps.size(); iii++)
	//		{
	//			if (reference[i].loci[ii].snps[iii] == 'n')
	//			{
	//				int this_locus = iii;
	//				reference[i].loci[ii].snps.erase(reference[i].loci[ii].snps.begin() + this_locus);
	//				reference[i].loci[ii].snp_freq.erase(reference[i].loci[ii].snp_freq.begin() + this_locus);
	//			}
	//		}
	//	}
	//	cout << "ref chrom " << reference[i].ID << " has " << reference[i].loci.size() << " loci.\n";
	//}

	ref_file.open(reference_out.str());
	ref_file << "Chrom\tLocus_name\tSNP_ID\tSNP_freq\n";
	for (i = 0; i < reference.size(); i++)
	{
		for (ii = 0; ii < reference[i].loci.size(); ii++)
		{
			ref_file << reference[i].ID << '\t' << reference[i].loci[ii].snp_name;
			for (iii = 0; iii < reference[i].loci[ii].snps.size(); iii++)
			{
				ref_file << '\t' << reference[i].loci[ii].snps[iii] << '\t' << reference[i].loci[ii].snp_freq[iii];
			}
			ref_file << '\n';
		}
	}
	ref_file.close();

	//set up necessary data structures to calculate allele frequencies
	for (i = 0; i < pops.size(); i++)
	{
		pop_afreqs_all.push_back(Genome(reference));
		exp_het_all.push_back(Genome(reference));
		obs_hets.push_back(Genome(reference));
	}

	for (v = 0; v < pops.size(); v++)
	{
		pops[v].establish_subsets(phen_of_interest, phen1_type, threshold);
		vector <Genome> temp;
		exp_het_sel.push_back(temp);
		obs_hets_s.push_back(temp);
		pop_afreqs_sel.push_back(temp);
		for (iv = 0; iv < pops[v].subsets.size(); iv++)
		{
			exp_het_sel[v].push_back(Genome(reference));
			obs_hets_s[v].push_back(Genome(reference));
			pop_afreqs_sel[v].push_back(Genome(reference));
		}
	}

	//within populations
	for (v = 0; v < pops.size(); v++)
	{
		cout << "\n\nCalculating allele frequencies for population " << pops[v].name << '\n';
		for (iv = 0; iv < pops[v].subsets.size(); iv++)
		{
			if (pops[v].adults.size() > 1)//allele frequency for the whole population
			{
				pop_afreqs_all[v] = pops[v].calc_allele_freq(pops[v].all_inds, reference);
				obs_hets[v] = pops[v].calc_obs_het(pops[v].all_inds, reference);
				exp_het_all[v] = pops[v].calc_exp_het(reference, pop_afreqs_all[v]);
			}
			if (pops[v].subsets[iv].size() > 1)//because the constructor creates 1 thing automatically
			{
				/*for (size_t test = 0; test < pops[v].subsets[iv].size(); test++)
				cout << pops[v].subsets[iv][test] << '\t';*/
				//cout << "there are " << pops[v].subsets[iv].size() << " individuals in subset " << iv << '\n';
				pop_afreqs_sel[v][iv] = pops[v].calc_allele_freq(pops[v].subsets[iv], reference);
				obs_hets_s[v][iv] = pops[v].calc_obs_het(pops[v].subsets[iv], reference);
				exp_het_sel[v][iv] = pops[v].calc_exp_het(reference, pop_afreqs_sel[v][iv]);
			}
		}
	}
	cout << "\nWriting within-population fst values to file.\n";
	summary_output.open(summary_name.str());
	summary_output << "Comparison\tNumGroup1\tNumGroup2\tSig99\tSig98\tSig95\tSig90\tSig80\tChr\tLocus";
	for (v = 0; v < pops.size(); v++)
	{
		for (iv = 0; iv < pops[v].subsets.size(); iv++)
		{
			for (iii = iv + 1; iii < pops[v].subsets.size(); iii++)
			{
				if (iv != iii)
				{
					//calc fsts between selected and non-selected, per population
					if (pops[v].subsets[iv].size() > 1 && pops[v].subsets[iii].size() > 1)
					{
						Genome fsts_in_pops(reference);
						Genome weighted_fsts_in_pops(reference);
						Cutoffs ci_cut_in_pops;
						fsts_in_pops = pops[v].Fst(pop_afreqs_sel[v][iv], pop_afreqs_sel[v][iii], reference);
						weighted_fsts_in_pops = pops[v].weighted_fst(fsts_in_pops, reference);
						ci_cut_in_pops = pops[v].calculate_fst_CIs(weighted_fsts_in_pops, reference);
						ofstream within_pop_fst;
						stringstream outfile_name;
						outfile_name << path_name << base_out_name << pop_names[v] << "_within_sub" << iv << "_sub" << iii << ".txt";
						within_pop_fst.open(outfile_name.str());
						within_pop_fst << "Chrom\tLocus\tGenDist\tObsHet\t";
						within_pop_fst << iii << "_ObsHet\t" << iv << "_ObsHet\tExpHet\t" << iii << "_ExpHet\t" << iv << "_ExpHet\t";
						within_pop_fst << "Fst\tWeightedFst\tSig99\tSig98\tSig95\tSig90\tSig80\t";
						within_pop_fst << iii << "MajAllele\t" << iii << "MajAllFreq\t" << iv << "MajAllele\t" << iv << "MajAllFreq\n";
						summary_output << '\n' << pops[v].name << "_sub" << iii << "-sub" << iv << '\t'
							<< pops[v].subsets[iii].size() - 1 << '\t' << pops[v].subsets[iv].size() - 1 << '\t';
						summary_output << ci_cut_in_pops.critical_99 << '\t' << ci_cut_in_pops.critical_98 << '\t'
							<< ci_cut_in_pops.critical_95 << '\t' << ci_cut_in_pops.critical_90 << '\t' << ci_cut_in_pops.critical_80;

						for (s_chr = 0; s_chr < reference.size(); s_chr++)
						{
							for (i = 0; i < reference[s_chr].loci.size(); i++)
							{
								within_pop_fst << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name << '\t'
									<< reference[s_chr].loci[i].distance << '\t';
								within_pop_fst << obs_hets[v].chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << obs_hets_s[v][iii].chrom_set[s_chr].loci[i].snp_freq[0] << '\t'
									<< obs_hets_s[v][iv].chrom_set[s_chr].loci[i].snp_freq[0] << '\t';
								within_pop_fst << exp_het_all[v].chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << exp_het_sel[v][iii].chrom_set[s_chr].loci[i].snp_freq[0] << '\t'
									<< exp_het_sel[v][iv].chrom_set[s_chr].loci[i].snp_freq[0] << '\t';
								within_pop_fst << fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0];
								within_pop_fst << '\t';
								if (weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_in_pops.critical_99)
								{
									within_pop_fst << 1 << '\t';
									if (stringency_value == 99)
										summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
								}
								else
									within_pop_fst << 0 << '\t';
								if (weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_in_pops.critical_98)
								{
									within_pop_fst << 1 << '\t';
									if (stringency_value == 98)
										summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
								}
								else
									within_pop_fst << 0 << '\t';
								if (weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_in_pops.critical_95)
								{
									within_pop_fst << 1 << '\t';
									if (stringency_value == 95)
										summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
								}
								else
									within_pop_fst << 0 << '\t';
								if (weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_in_pops.critical_90)
								{
									within_pop_fst << 1 << '\t';
									if (stringency_value == 90)
										summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
								}
								else
									within_pop_fst << 0 << '\t';
								if (weighted_fsts_in_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_in_pops.critical_80)
								{
									within_pop_fst << 1;
									if (stringency_value == 80)
										summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
								}
								else
									within_pop_fst << 0;
								maj_all_freq_1 = 0;
								maj_all_freq_2 = 0;
								for (ii = 0; ii < reference[s_chr].loci[i].snps.size(); ii++)
								{
									if (pop_afreqs_sel[v][iii].chrom_set[s_chr].loci[i].snp_freq[ii] > maj_all_freq_1)
									{
										maj_all_freq_1 = pop_afreqs_sel[v][iii].chrom_set[s_chr].loci[i].snp_freq[ii];
										maj_all_loc_1 = ii;
									}
									if (pop_afreqs_sel[v][iv].chrom_set[s_chr].loci[i].snp_freq[ii] > maj_all_freq_2)
									{
										maj_all_freq_2 = pop_afreqs_sel[v][iv].chrom_set[s_chr].loci[i].snp_freq[ii];
										maj_all_loc_2 = ii;
									}
								}
								within_pop_fst << '\t' << pop_afreqs_sel[v][iii].chrom_set[s_chr].loci[i].snps[maj_all_loc_1] << '\t'
									<< maj_all_freq_1 << '\t';
								within_pop_fst << pop_afreqs_sel[v][iv].chrom_set[s_chr].loci[i].snps[maj_all_loc_2] << '\t'
									<< maj_all_freq_2;
								within_pop_fst << '\n';
							}
						}
						within_pop_fst.close();
						cout << "Population " << pop_names[v] << " fst data written to file.\n\n";
					}
				}
			}
		}
	}

	//compare populations
	int compare_count = 0;
	for (v = 0; v < pops.size(); v++)
	{
		pop_afreqs_all[v] = pops[v].calc_allele_freq(pops[v].all_inds, reference);
		obs_hets[v] = pops[v].calc_obs_het(pops[v].all_inds, reference);
		exp_het_all[v] = pops[v].calc_exp_het(reference, pop_afreqs_all[v]);
	}

	for (v = 0; v < pops.size(); v++)
	{
		for (iv = v + 1; iv < pops.size(); iv++)
		{
			if (v != iv)
			{
				if (pops[v].adults.size() > 0 && pops[iv].adults.size() > 0)
				{
					ofstream between_pop_fst;
					out_name = path_name + base_out_name + pop_names[v] + "_vs_" + pop_names[iv] + ".txt";
					between_pop_fst.open(out_name);
					between_pop_fst << "Chrom\tLocus\tDist\tObsHetPop1\tObsHetPop2\tExpHetPop1\tExpHetPop2\tFst\tWeightedFst\t"
						<< "Sig99\tSig98\tSig95\tSig90\tSig80\tMajAllele1\tMajAllFreq1\tMajAllele2\tMajAllFreq2\n";

					cout << "\nCalculating Fst values between " << pops[v].name << " and " << pops[iv].name << ".\n";
					Genome fsts_btwn_pops(reference);
					Genome weighted_fsts_btwn_pops(reference);
					Cutoffs ci_cut_btwn_pops;

					fsts_btwn_pops = pops[v].Fst(pop_afreqs_all[v], pop_afreqs_all[iv], reference);
					weighted_fsts_btwn_pops = pops[v].weighted_fst(fsts_btwn_pops, reference);
					ci_cut_btwn_pops = pops[v].calculate_fst_CIs(weighted_fsts_btwn_pops, reference);
					compare_count++;

					summary_output << '\n' << pops[v].name << "-" << pops[iv].name << '\t'
						<< pops[v].adults.size() - 1 << '\t' << pops[iv].adults.size() - 1 << '\t';
					summary_output << ci_cut_btwn_pops.critical_99 << '\t' << ci_cut_btwn_pops.critical_98 << '\t'
						<< ci_cut_btwn_pops.critical_95 << '\t' << ci_cut_btwn_pops.critical_90 << '\t' << ci_cut_btwn_pops.critical_80;


					//write to file
					for (s_chr = 0; s_chr < reference.size(); s_chr++)
					{
						for (i = 0; i < reference[s_chr].loci.size(); i++)
						{
							between_pop_fst << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name << '\t'
								<< reference[s_chr].loci[i].distance << '\t';
							between_pop_fst << obs_hets[v].chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << obs_hets[iv].chrom_set[s_chr].loci[i].snp_freq[0] << '\t';
							between_pop_fst << exp_het_all[v].chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << exp_het_all[iv].chrom_set[s_chr].loci[i].snp_freq[0] << '\t';
							between_pop_fst << fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] << '\t' << weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0];
							if (weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_btwn_pops.critical_99)
							{
								between_pop_fst << '\t' << 1;
								if (stringency_value == 99)
									summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
							}
							else
								between_pop_fst << '\t' << 0;
							if (weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_btwn_pops.critical_98)
							{
								between_pop_fst << '\t' << 1;
								if (stringency_value == 98)
									summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
							}
							else
								between_pop_fst << '\t' << 0;
							if (weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_btwn_pops.critical_95)
							{
								between_pop_fst << '\t' << 1;
								if (stringency_value == 95)
									summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
							}
							else
								between_pop_fst << '\t' << 0;
							if (weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_btwn_pops.critical_90)
							{
								between_pop_fst << '\t' << 1;
								if (stringency_value == 90)
									summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
							}
							else
								between_pop_fst << '\t' << 0;
							if (weighted_fsts_btwn_pops.chrom_set[s_chr].loci[i].snp_freq[0] >= ci_cut_btwn_pops.critical_80)
							{
								between_pop_fst << '\t' << 1;
								if (stringency_value == 80)
									summary_output << '\t' << reference[s_chr].ID << '\t' << reference[s_chr].loci[i].snp_name;
							}
							else
								between_pop_fst << '\t' << 0;
							maj_all_freq_1 = 0;
							maj_all_freq_2 = 0;
							for (ii = 0; ii < reference[s_chr].loci[i].snps.size(); ii++)
							{
								if (pop_afreqs_all[v].chrom_set[s_chr].loci[i].snp_freq[ii] > maj_all_freq_1)
								{
									maj_all_freq_1 = pop_afreqs_all[v].chrom_set[s_chr].loci[i].snp_freq[ii];
									maj_all_loc_1 = ii;
								}
								if (pop_afreqs_all[iv].chrom_set[s_chr].loci[i].snp_freq[ii] > maj_all_freq_2)
								{
									maj_all_freq_2 = pop_afreqs_all[iv].chrom_set[s_chr].loci[i].snp_freq[ii];
									maj_all_loc_2 = ii;
								}
							}
							between_pop_fst << '\t' << pop_afreqs_all[v].chrom_set[s_chr].loci[i].snps[maj_all_loc_1] << '\t'
								<< maj_all_freq_1 << '\t';
							between_pop_fst << pop_afreqs_all[iv].chrom_set[s_chr].loci[i].snps[maj_all_loc_2] << '\t'
								<< maj_all_freq_2;
							between_pop_fst << '\n';
						}
					}
				}
			}
		}
	}
	cout << "\nPerformed " << compare_count << " pairwise population-level Fst comparisons.\n";
	summary_output.close();
	allele_out.open(allele_name.str());
	allele_out << "Population";
	for (ii = 0; ii < reference.size(); ii++)
	{
		for (iii = 0; iii < reference[ii].loci.size(); iii++)
		{
			allele_out << '\t' << ii << "." << iii << "_allele\t" << ii << "." << iii << "_freq";
		}
	}
	char maj_allele;
	for (i = 0; i < pops.size(); i++)
	{
		allele_out << '\n' << pops[i].name;
		for (ii = 0; ii < reference.size(); ii++)
		{
			for (iii = 0; iii < reference[ii].loci.size(); iii++)
			{
				maj_all_freq_1 = 0;
				for (iv = 0; iv < reference[ii].loci[iii].snps.size(); iv++)
				{
					if (pop_afreqs_all[i].chrom_set[ii].loci[iii].snp_freq[iv] > maj_all_freq_1)
					{
						maj_all_freq_1 = pop_afreqs_all[i].chrom_set[ii].loci[iii].snp_freq[iv];
						maj_allele = pop_afreqs_all[i].chrom_set[ii].loci[iii].snps[iv];
					}
				}
				allele_out << '\t' << maj_allele << '\t' << maj_all_freq_1;
			}
		}
	}
	allele_out.close();
	//end of the program
	if (interactivemode)
	{
		cout << "\nDone! Enter integer to quit.\n";
		cin >> end;
	}

	return 0;
}