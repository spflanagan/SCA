//Date: 1 February 2015
//Author: Sarah P. Flanagan
//Purpose of this program: to generate a dataset with known allele frequencies to test my sparrow analysis program
//This will output an individual information file, a .ped file, a .map file, and an allele frequencies file

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "rand_num_distributions.h"
#include "classes.h"

using namespace std;

int main()
{
	int end, i, ii, iii;
	size_t j, jj, jjj, iv, v;
	int num_chrom, max_alleles_per_locus, max_loci_per_chrom, num_individuals, max_types_phen1, num_populations, max_phen2;
	vector <char> all_snps = { 'A', 'C', 'T', 'G' };
	vector <char> available_snps;
	vector <Genome> allele_freqs;
	vector <Chromosome> reference;
	vector <string> ind_ids;
	char temp_snp;
	int temp_location, snp_count, location_counter, num_alleles, rand_num, rand1, rand2, num_years;
	ofstream map_file, ped_file, ind_file, allele_file;
	stringstream map_name, ped_name, ind_name, allele_name;
	string base_out_name;
	double average_loci_per_chrom;

	num_chrom = 10;
	max_loci_per_chrom = 100;
	average_loci_per_chrom = 0.75;//most will be about 0.75*max_loci_per_chrom
	max_alleles_per_locus = 3;
	num_individuals = 300;
	max_types_phen1 = 4;
	num_populations = 3;
	max_phen2 = 20;
	num_years = 5;
	base_out_name = "test_data";
	map_name << base_out_name << ".map";
	ped_name << base_out_name << ".ped";
	ind_name << base_out_name << ".ind.txt";
	allele_name << base_out_name << ".alleles.txt";

	cout << "Generating .map file\n";
	snp_count = 0;
	for (i = 0; i < num_chrom; i++)
	{
		//set up reference
		reference.push_back(Chromosome());
		reference[i].ID = i;
		location_counter = 0;
		for (ii = 0; ii < max_loci_per_chrom; ii++)
		{
			if (genrand() < average_loci_per_chrom)
			{
				reference[i].loci.push_back(Locus_info());
				stringstream locus_id;
				locus_id << num_chrom << location_counter;
				reference[i].loci[location_counter].snp_name = locus_id.str();
				available_snps = { 'A', 'C', 'T', 'G' };
				num_alleles = randnum(max_alleles_per_locus) + 1;
				for (iii = 0; iii < num_alleles; iii++)
				{
					temp_location = randnum(available_snps.size());
					temp_snp = available_snps[temp_location];
					reference[i].loci[location_counter].snps.push_back(temp_snp);
					reference[i].loci[location_counter].snp_freq.push_back(0);
					snp_count++;
					available_snps.erase(available_snps.begin() + temp_location);
				}
				location_counter++;
			}
		}
	}//end reference setup
	for (j = 0; j < reference.size(); j++)
	{
		reference[j].remove_defaults();
		cout << "Chromsome " << reference[j].ID << " has " << reference[j].loci.size() << " loci.\n";
	}

	map_file.open(map_name.str());
	//map file has four columns: chromosome, SNP ID, genetic distance, bp position of snp
	for (j = 0; j < reference.size(); j++)
	{
		for (jj = 0; jj < reference[j].loci.size(); jj++)
		{
			map_file << reference[j].ID << '\t' << reference[j].loci[jj].snp_name << '\t' 
				<< j + jj << '\t' << jj << '\n';
		}
	}
	map_file.close();


	//create individuals
	vector <Population> pops;
	for (i = 0; i < num_populations; i++)
	{
		pops.push_back(Population());
		stringstream name;
		name << "pop" << i;
		pops[i].name = name.str();
	}
	for (i = 0; i < num_individuals; i++)
	{
		stringstream name;
		name << i + 1;
		ind_ids.push_back(name.str());
	}
	location_counter = 0;
	for (i = 0; i < num_individuals; i++)
	{
		rand_num = randnum(num_populations);
		for (ii = 0; ii < num_populations; ii++)
		{
			if (rand_num == ii)
			{
				pops[ii].adults.push_back(Individual());
				pops[ii].adults[pops[ii].adults.size() - 1].ID = ind_ids[location_counter];
				location_counter++;
				if (genrand() < 0.5)
					pops[ii].adults[pops[ii].adults.size() - 1].female = true;
				else
					pops[ii].adults[pops[ii].adults.size() - 1].female = false;
				rand2 = randnum(num_years);
				for (iii = 1; iii < rand2; iii++)//1 is already initialized in Individual()
				{
					pops[ii].adults[pops[ii].adults.size() - 1].year.push_back(iii);
					pops[ii].adults[pops[ii].adults.size() - 1].phenotype1.push_back(randnum(max_types_phen1));
					pops[ii].adults[pops[ii].adults.size() - 1].phenotype2.push_back(randnum(max_phen2));
					if (genrand() < 0.5 && pops[ii].adults[pops[ii].adults.size() - 1].phenotype3[iii-1] == true)
						pops[ii].adults[pops[ii].adults.size() - 1].phenotype3.push_back(true);
					else
						pops[ii].adults[pops[ii].adults.size() - 1].phenotype3.push_back(false);
				}
				//assign genotypes
				for (jj = 0; jj < reference.size(); jj++)
				{
					pops[ii].adults[pops[ii].adults.size() - 1].allele1.push_back(Chromosome());
					pops[ii].adults[pops[ii].adults.size() - 1].allele2.push_back(Chromosome());
					for (jjj = 0; jjj < reference[jj].loci.size(); jjj++)
					{
						rand1 = randnum(reference[jj].loci[jjj].snps.size()); 
						rand2 = randnum(reference[jj].loci[jjj].snps.size());
						pops[ii].adults[pops[ii].adults.size() - 1].allele1[jj].loci.push_back(Locus_info());
						pops[ii].adults[pops[ii].adults.size() - 1].allele2[jj].loci.push_back(Locus_info());
						pops[ii].adults[pops[ii].adults.size() - 1].allele1[jj].loci[jjj].snps[0] = reference[jj].loci[jjj].snps[rand1];
						pops[ii].adults[pops[ii].adults.size() - 1].allele2[jj].loci[jjj].snps[0] = reference[jj].loci[jjj].snps[rand2];
						pops[ii].adults[pops[ii].adults.size() - 1].allele1[jj].loci[jjj].snp_freq[0] = 0;
						pops[ii].adults[pops[ii].adults.size() - 1].allele2[jj].loci[jjj].snp_freq[0] = 0;
						pops[ii].adults[pops[ii].adults.size() - 1].allele1[jj].loci[jjj].snp_name = reference[jj].loci[jjj].snp_name;
						pops[ii].adults[pops[ii].adults.size() - 1].allele2[jj].loci[jjj].snp_name = reference[jj].loci[jjj].snp_name;
					}
				}
				for (jj = 0; jj < reference.size(); jj++)
				{
					pops[ii].adults[pops[ii].adults.size() - 1].allele1[jj].remove_defaults();
					pops[ii].adults[pops[ii].adults.size() - 1].allele2[jj].remove_defaults();
				}
			}
		}
	}
	for (i = 0; i < num_populations; i++)
		pops[i].remove_defaults();

	//generate ped file: Fam ID, Ind ID, Pat ID, Mat ID, sex, phenotype (not a real one)
	cout << "\nWriting to .ped file.\n";
	ped_file.open(ped_name.str());
	for (i = 0; i < num_populations; i++)
	{
		for (j = 0; j < pops[i].adults.size(); j++)
		{
			ped_file << 1 << '\t' << pops[i].adults[j].ID << '\t' << 0 << '\t' << 0 << '\t';
			if (pops[i].adults[j].female)
				ped_file << 0 << '\t' << -9;
			if (!pops[i].adults[j].female)
				ped_file << 1 << '\t' << -9;
			for (jj = 0; jj < reference.size(); jj++)
			{
				for (jjj = 0; jjj < reference[jj].loci.size(); jjj++)
				{
					ped_file << '\t' << pops[i].adults[j].allele1[jj].loci[jjj].snps[0]
						<< '\t' << pops[i].adults[j].allele2[jj].loci[jjj].snps[0];
				}
			}
			ped_file << '\n';
		}
	}
	ped_file.close();

	//write individual info file
	cout << "\nWriting .ind.txt file\n";
	ind_file.open(ind_name.str());
	ind_file << "ID\tYear\tLocality\tPop\tSex\tAge\tPhen1\tPhen2\tPhen3";
	int count = 0;
	for (i = 0; i < num_populations; i++)
	{
		for (j = 0; j < pops[i].adults.size(); j++)
		{
			for (jj = 0; jj < pops[i].adults[j].year.size(); jj++)
			{
				ind_file << '\n';
				ind_file << pops[i].adults[j].ID << '\t' << jj + 1 << '\t' << 0 << '\t' << pops[i].name << '\t';
				if (pops[i].adults[j].female)
					ind_file << 0;
				if (!pops[i].adults[j].female)
					ind_file << 1;
				ind_file << '\t' << jj + 2 << '\t' << pops[i].adults[j].phenotype1[jj] << '\t' << pops[i].adults[j].phenotype2[jj]
					<< '\t' << pops[i].adults[j].phenotype3[jj];
				count++;
			}
		}
	}
	ind_file.close();
	cout << '\n';
	for (i = 0; i < num_populations; i++)
	{
		cout << "Population " << pops[i].name << " has " << pops[i].adults.size() << " individuals.\n";
		allele_freqs.push_back(Genome(reference));
		for (j = 0; j < reference.size(); j++)
		{
			//allele_freqs[i].chrom_set.push_back(Chromosome());
			for (jj = 0; jj < reference[j].loci.size(); jj++)
			{
				//allele_freqs[i].chrom_set[j].loci.push_back(Locus_info());
				allele_freqs[i].chrom_set[j].ID = reference[j].ID;
				allele_freqs[i].chrom_set[j].loci[jj].snp_name = reference[j].loci[jj].snp_name;
				for (jjj = 0; jjj < reference[j].loci[jj].snps.size(); jjj++)
				{
					allele_freqs[i].chrom_set[j].loci[jj].snps.push_back(reference[j].loci[jj].snps[jjj]);
					allele_freqs[i].chrom_set[j].loci[jj].snp_freq.push_back(0);
				}
			}
			allele_freqs[i].chrom_set[j].remove_defaults();
		}
	}

	for (i = 0; i < num_populations; i++)
	{
		for (j = 0; j < pops[i].adults.size(); j++)
		{
			for (jj = 0; jj < reference.size(); jj++)
			{
				for (jjj = 0; jjj < reference[jj].loci.size(); jjj++)
				{
					for (iv = 0; iv < reference[jj].loci[jjj].snps.size(); iv++)
					{
						if (pops[i].adults[j].allele1[jj].loci[jjj].snps[0] == reference[jj].loci[jjj].snps[iv])
							allele_freqs[i].chrom_set[jj].loci[jjj].snp_freq[iv] = allele_freqs[i].chrom_set[jj].loci[jjj].snp_freq[iv] + 1;
						if (pops[i].adults[j].allele2[jj].loci[jjj].snps[0] == reference[jj].loci[jjj].snps[iv])
							allele_freqs[i].chrom_set[jj].loci[jjj].snp_freq[iv] = allele_freqs[i].chrom_set[jj].loci[jjj].snp_freq[iv] + 1;
					}
					
				}
			}
		}
	}
	cout << "\nWriting allele frequencies to file.\n";
	allele_file.open(allele_name.str());
	allele_file << "Population";
	for (j = 0; j < reference.size(); j++)
	{
		for (jj = 0; jj < reference[j].loci.size(); jj++)
		{
			allele_file << '\t' << j << "." << jj << "_allele" << '\t' << j << "." << jj << "_freq";
		}
	}
	double maj_allele_freq;
	char maj_allele;
	for (i = 0; i < num_populations; i++)
	{
		allele_file << '\n' << pops[i].name;
		for (j = 0; j < allele_freqs[i].chrom_set.size(); j++)
		{
			for (jj = 0; jj < allele_freqs[i].chrom_set[j].loci.size(); jj++)
			{
				for (jjj = 0; jjj < allele_freqs[i].chrom_set[j].loci[jj].snps.size(); jjj++)
				{
					double adult_num = pops[i].adults.size();
					allele_freqs[i].chrom_set[j].loci[jj].snp_freq[jjj] = allele_freqs[i].chrom_set[j].loci[jj].snp_freq[jjj] /
						(2 * adult_num);
				}
				maj_allele_freq = 0;
				for (jjj = 0; jjj < allele_freqs[i].chrom_set[j].loci[jj].snps.size(); jjj++)
				{
					if (allele_freqs[i].chrom_set[j].loci[jj].snp_freq[jjj] > maj_allele_freq)
					{
						maj_allele_freq = allele_freqs[i].chrom_set[j].loci[jj].snp_freq[jjj];
						maj_allele = allele_freqs[i].chrom_set[j].loci[jj].snps[jjj];
					}
				}
				allele_file << '\t' << maj_allele << '\t' << maj_allele_freq;
			}
		}
	}
	allele_file.close();
	cout << "\nDone! Enter integer to quit.\n";
	cin >> end;
	return 0;
}