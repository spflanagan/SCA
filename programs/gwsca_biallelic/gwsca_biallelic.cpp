//Author: Sarah P. Flanagan
//Date: 26 January 2016
//Purpose: Take reformatted stacks output from many individuals, match it with
//information about their life history stage, and generate Fst values between many 
//different comparisons
//e.g., male-female, adult-offspring, inferred mothers-collected females
//This uses biallelic SNP data gleaned from stacks snps and matches files.
//This program is basically the same as gwsca_haplotypes.cpp but it the input files don't have counts
//The ind_info file it reads has four columns, with one row for each individual:
//1. file names (with paths to the location)
//2. The individual ID
//3. Sex
//4. Status (e.g. pregnant/nonpregnant, mated/nonmated, adult/offspring)

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
	string filename, ind_id, sex, age, status;

	individual()
	{
		filename = ind_id = sex = age = status = string();
	}
};

class locus
{
public:
	string snp_id;
	int count, pos, col, cat_id;
	vector<string> alleles;

	locus()
	{
		count = int();
		pos = int();
		col = int();
		cat_id = int();
		snp_id = string();
		alleles = vector < string >();
	}
};

class locus_statistics
{
public:
	bool haploid;
	string group_name;
	vector<vector<double>> allele_frequency;
	vector<double> hs, ho;
	vector<int> major_allele_index;
	int two_n;

	locus_statistics()
	{
		haploid = bool();
		allele_frequency = vector<vector<double>>();
		hs = vector<double>();
		ho = vector<double>();
		group_name = string();
		two_n = int();
		major_allele_index = vector<int>();
	}
};

int main()
{
	int i, ii, iii, end, pop_index, al1_count, al2_count, loc_id, index, al1i, al2i, removed_count, bp, col, cat_id;
	size_t t, tt, ttt;
	string line, filename, ind_id, snp_id, sex, age, status, allele1, allele2, chr;
	string ind_info_name, fst_out_name, summary_out_name, alleles_out_name, sumstats_name;
	ifstream ind_info, sumstats;
	ofstream fst_out, summary_out, alleles_out;
	vector<individual> inds;
	vector<locus> reference;
	vector<locus_statistics> pop_stats;

	sumstats_name = "../../results/stacks/batch_1.sumstats.tsv";
	ind_info_name = "../../results/biallelic/ind_info_gwsca.txt";
	ind_info.open(ind_info_name);
	FileTest(ind_info, ind_info_name);
	while (universal_getline(ind_info, line))
	{
		if (!ind_info.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> filename >> ind_id >> sex >> age >> status;
			inds.push_back(individual());
			inds.back().filename = filename;
			inds.back().ind_id = ind_id;
			inds.back().age = age;
			inds.back().status = status;
			inds.back().sex = sex;
		}
	}
	ind_info.close();

	cout << "\nRead in data for " << inds.size() << " individuals.\n";

	//establish the population statistics counters
	bool age_found, sex_found, status_found;
	for (t = 0; t < inds.size(); t++)
	{
		age_found = sex_found = status_found = false;
		for (tt = 0; tt < pop_stats.size(); tt++)
		{
			if (inds[t].sex == pop_stats[tt].group_name)
			{
				sex_found = true;
				pop_stats[tt].two_n++;
			}
			if (inds[t].age == pop_stats[tt].group_name)
			{
				age_found = true;
				pop_stats[tt].two_n++;
			}
			if (inds[t].status == pop_stats[tt].group_name)
			{
				status_found = true;
				pop_stats[tt].two_n++;
			}
		}
		if (age_found == false)
		{
			if (inds[t].age != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].age;
				pop_stats.back().two_n = 1;
			}
		}
		if (sex_found == false)
		{
			if (inds[t].sex != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].sex;
				pop_stats.back().two_n = 1;
			}
		}
		if (status_found == false)
		{
			if (inds[t].status != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].status;
				pop_stats.back().two_n = 1;
			}
		}
	}
	cout << "\nEstablished " << pop_stats.size() << " groups:";
	for (t = 0; t < pop_stats.size(); t++)
	{
		cout << '\n' << pop_stats[t].group_name;
	}
	cout << "\n\n";

	//let's make a reference, using sumstats.tsv
	sumstats.open(sumstats_name);
	FileTest(sumstats, sumstats_name);
	int line_count = 0;
	while (universal_getline(sumstats, line))
	{
		if (!sumstats.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> i >> cat_id >> chr >> bp >> col;
				reference.push_back(locus());
				reference.back().cat_id = cat_id;
				reference.back().col = col;
				reference.back().pos = bp;
				stringstream sstemp;
				sstemp << cat_id << "." << col;
				reference.back().snp_id = sstemp.str();
				reference.back().count = 0;
				line_count++;
				//also start the allele frequency counters
				for (ttt = 0; ttt < pop_stats.size(); ttt++)
				{
					pop_stats[ttt].major_allele_index.push_back(int());
					pop_stats[ttt].allele_frequency.push_back(vector<double>());
					pop_stats[ttt].allele_frequency.back().push_back(0);
					pop_stats[ttt].allele_frequency.back().push_back(0);
					pop_stats[ttt].ho.push_back(0);
					pop_stats[ttt].hs.push_back(1);
					pop_stats[ttt].haploid = false;
				}
			}
		}
	}
	sumstats.close();
	cout << "\nThe reference contains " << line_count << " SNPs.\n";

	//now need to read in the allele data and add alleles to reference
	for (t = 0; t < inds.size(); t++)
	{
		vector<int> indices;
		for (tt = 0; tt < pop_stats.size(); tt++)
		{
			if (inds[t].sex == pop_stats[tt].group_name)
				indices.push_back(tt);
			if (inds[t].age == pop_stats[tt].group_name)
				indices.push_back(tt);
			if (inds[t].status == pop_stats[tt].group_name)
				indices.push_back(tt);
		}
		ifstream ind_file;
		ind_file.open(inds[t].filename);
		FileTest(ind_file, inds[t].filename);
		int line_count = 0;
		removed_count = 0;
		while (universal_getline(ind_file, line))
		{
			if (!ind_file.eof())
			{
				if (line_count > 0)
				{
					stringstream ss;
					ss.str(line);
					ss >> snp_id >> allele1 >>  allele2;
					index = -5;
					for (tt = 0; tt < reference.size(); tt++)
					{
						if (inds[t].sex == "MOM")
						{
							string stemp;
							stemp = find_and_replace(snp_id, ".", "\t");
							stringstream momss;
							momss.str(stemp);
							int cat, pos;
							momss >> cat >> pos;

							if (reference[tt].cat_id == cat && reference[tt].pos)
								index = tt;
						}
						else
						{
							if (reference[tt].snp_id == snp_id)
								index = tt;
						}
					}
					if (index >= 0)
					{//the locus is in the reference
						reference[index].count++;
						al1i = al2i = -5;
						for (ttt = 0; ttt < reference[index].alleles.size(); ttt++)
						{//is the allele in the reference?
							if (reference[index].alleles[ttt] == allele1)
								al1i = ttt;
							if (reference[index].alleles[ttt] == allele2)
								al2i = ttt;
						}
						if (al1i < 0 && allele1 != "-")
						{
							reference[index].alleles.push_back(allele1);
							
								al1i = reference[index].alleles.size() - 1;
						}
						if (al2i < 0 && allele2 != "-")
						{
							if (allele2 != "-" && allele1 != allele2)
							{
								reference[index].alleles.push_back(allele2);
								al2i = reference[index].alleles.size() - 1;
							}
							if (allele1 == allele2)//it's homozygous
								al2i = al1i;
						}
					
						//check if it's more than two alleles
						if (reference[index].alleles.size() > 2)
						{//remove it from the records
							cout << "\nAllele " << index << " had more than two alleles. Removing it from dataset.";
							reference.erase(reference.begin() + index);
							for (i = 0; i < pop_stats.size(); i++)
								pop_stats[i].allele_frequency.erase(pop_stats[i].allele_frequency.begin() + index);
							removed_count++;
						}
						else
						{	//add the counts to the index
							for (i = 0; i < pop_stats.size(); i++)
								pop_stats[i].allele_frequency[index].push_back(0);
							for (ttt = 0; ttt < indices.size(); ttt++)
							{
								pop_stats[indices[ttt]].allele_frequency[index][al1i]++;
								if (allele2 != "-")
									pop_stats[indices[ttt]].allele_frequency[index][al2i]++;
								if (allele1 == allele2)
									pop_stats[indices[ttt]].ho[index]++;
								if (allele2 == "-")
									pop_stats[indices[ttt]].haploid = true;
							}
						}
					}
					else
					{// the locus is not in the reference
						//cout << "Locus " << snp_id << " is not in the reference.\n";
					}
				}
				line_count++;
			}
		}
		ind_file.close();
	}
	
	cout << "\nRemoved " << removed_count << " SNP records that were not biallelic.\n";
	cout << "Recorded " << reference.size() << " biallelic SNPs in " << inds.size() << " individuals spread among " << pop_stats.size() << " groups.\n";
	summary_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_summary.txt";
	summary_out.open(summary_out_name);
	summary_out << "Pop\tLocusName\tHs\tHo\tNumAlleles";
	vector<vector<double>> fsts;
	vector<vector<double>> ht;
	for (t = 0; t < pop_stats.size(); t++)
	{
		for (tt = t + 1; tt < pop_stats.size(); tt++)
		{
			if (t != tt)
			{
				fsts.push_back(vector<double>());
				ht.push_back(vector<double>());
				for (i = 0; i < pop_stats[t].allele_frequency.size(); i++)
				{
					fsts.back().push_back(0);
					ht.back().push_back(1);
				}
			}
		}
	}
	//now calculate frequencies
	for (t = 0; t < pop_stats.size(); t++)
	{
		//cout << "\n\n" << pop_stats[t].group_name;
		if (pop_stats[t].two_n > 0)
		{
			if (!pop_stats[t].haploid)
				pop_stats[t].two_n = pop_stats[t].two_n * 2;
			for (tt = 0; tt < pop_stats[t].allele_frequency.size(); tt++)
			{
				pop_stats[t].hs[tt] = 1;
				if (pop_stats[t].haploid)
					pop_stats[t].ho[tt] = pop_stats[t].ho[tt] / pop_stats[t].two_n;
				else
					pop_stats[t].ho[tt] = pop_stats[t].ho[tt] / (pop_stats[t].two_n / 2);
				for (ttt = 0; ttt < pop_stats[t].allele_frequency[tt].size(); ttt++)
				{
					pop_stats[t].allele_frequency[tt][ttt] = pop_stats[t].allele_frequency[tt][ttt] / pop_stats[t].two_n;
					pop_stats[t].hs[tt] = pop_stats[t].hs[tt] - (pop_stats[t].allele_frequency[tt][ttt] * pop_stats[t].allele_frequency[tt][ttt]);
				}

				summary_out << '\n' << pop_stats[t].group_name << '\t' << tt << '\t' << pop_stats[t].hs[tt] << '\t' << pop_stats[t].ho[tt] << '\t' << pop_stats[t].allele_frequency[tt].size();
			}
		}
		else
			cout << pop_stats[t].group_name << " has no members.\n";
	}
	summary_out.close();

	cout << "\nWriting Allele Information to File\n";
	alleles_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_haplotypes_alleles.txt";
	alleles_out.open(alleles_out_name);
	alleles_out << "Locus\tAllele";
	for (t = 0; t < pop_stats.size(); t++)
		alleles_out << "\t" << pop_stats[t].group_name << "Count";
	for (t = 0; t < reference.size(); t++)
	{
		for (tt = 0; tt < reference[t].alleles.size(); tt++)
		{
			alleles_out << '\n' << reference[t].cat_id << '\t' << reference[t].alleles[tt];
			for (ttt = 0; ttt < pop_stats.size(); ttt++)
				alleles_out << '\t' << pop_stats[ttt].allele_frequency[t][tt] * pop_stats[ttt].two_n;
		}
	}
	alleles_out.close();

	index = 0;
	for (t = 0; t < pop_stats.size(); t++)
	{
		for (tt = t + 1; tt < pop_stats.size(); tt++)
		{
			if (t != tt)
			{
				for (i = 0; i < reference.size(); i++)
				{
					ht[index][i] = 1;
					double avgp = 0;
					for (ii = 0; ii < pop_stats[t].allele_frequency[i].size(); ii++)
					{
						avgp = (pop_stats[t].allele_frequency[i][ii] + pop_stats[tt].allele_frequency[i][ii]) / 2;
						ht[index][i] = ht[index][i] - (avgp * avgp);
					}
				}
				index++;
			}
		}
	}
	index = 0;

	cout << "\nWriting Fst Values to File.\n";
	fst_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_fsts.txt";
	fst_out.open(fst_out_name);
	fst_out << "Locus";
	for (t = 0; t < pop_stats.size(); t++)
	{
		for (tt = t + 1; tt < pop_stats.size(); tt++)
		{
			fst_out << '\t' << pop_stats[t].group_name << "-" << pop_stats[tt].group_name;
		}
	}
	for (i = 0; i < reference.size(); i++)
	{
		fst_out << '\n' << reference[i].cat_id;
		index = 0;
		for (t = 0; t < pop_stats.size(); t++)
		{
			for (tt = t + 1; tt < pop_stats.size(); tt++)
			{
				if (t != tt)
				{
					double hs = (pop_stats[t].hs[i] + pop_stats[tt].hs[i]) / 2;
					if (ht[index][i] > 0)
						fsts[index][i] = (ht[index][i] - hs) / ht[index][i];
					else
						fsts[index][i] = -1;
					fst_out << '\t' << fsts[index][i];
					index++;
				}
			}
		}
	}
	fst_out.close();

	//calculate lumped Fsts
	for (t = 0; t < pop_stats.size(); t++)
	{
		for (tt = 0; tt < reference.size(); tt++)
		{
			double max_freq = 0;
			for (ttt = 0; ttt < pop_stats[t].allele_frequency[tt].size(); ttt++)
			{
				if (pop_stats[t].allele_frequency[tt][ttt] > max_freq)
				{
					max_freq = pop_stats[t].allele_frequency[tt][ttt];
					pop_stats[t].major_allele_index[tt] = ttt;
				}
			}
		}
	}

	cout << "\nWriting Lumped Fst Values to File.\n";
	fst_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_fsts_lumped.txt";
	fst_out.open(fst_out_name);
	fst_out << "Locus";
	for (t = 0; t < pop_stats.size(); t++)
	{
		for (tt = t + 1; tt < pop_stats.size(); tt++)
		{
			fst_out << '\t' << pop_stats[t].group_name << "-" << pop_stats[tt].group_name;
		}
	}
	double hs1, hs2, ht_temp, fst, af_sum1, af_sum2, avgm, avgo;
	for (i = 0; i < reference.size(); i++)
	{
		fst_out << '\n' << reference[i].cat_id;
		for (t = 0; t < pop_stats.size(); t++)
		{
			for (tt = t + 1; tt < pop_stats.size(); tt++)
			{
				if (t != tt)
				{
					af_sum1 = 0;
					for (ttt = 0; ttt < pop_stats[t].allele_frequency[i].size(); ttt++)
					{
						if (ttt != pop_stats[t].major_allele_index[tt])
							af_sum1 = af_sum1 + pop_stats[t].allele_frequency[i][ttt];
					}
					hs1 = 1 - ((af_sum1*af_sum1) +
						(pop_stats[t].allele_frequency[i][pop_stats[t].major_allele_index[i]] * pop_stats[t].allele_frequency[i][pop_stats[t].major_allele_index[i]]));

					af_sum2 = 0;
					for (ttt = 0; ttt < pop_stats[tt].allele_frequency[i].size(); ttt++)
					{
						if (ttt != pop_stats[tt].major_allele_index[tt])
							af_sum2 = af_sum2 + pop_stats[tt].allele_frequency[i][ttt];
					}
					hs2 = 1 - ((af_sum2*af_sum2) +
						(pop_stats[tt].allele_frequency[i][pop_stats[tt].major_allele_index[i]] * pop_stats[tt].allele_frequency[i][pop_stats[tt].major_allele_index[i]]));

					avgm = (pop_stats[t].allele_frequency[i][pop_stats[tt].major_allele_index[i]] + pop_stats[tt].allele_frequency[i][pop_stats[tt].major_allele_index[i]]) / 2;
					avgo = (af_sum1 + af_sum2) / 2;
					ht_temp = 1 - ((avgm * avgm) + (avgo*avgo));
					double hs = (hs1 + hs2) / 2;
					if (ht_temp > 0)
						fst = (ht_temp - hs) / ht_temp;
					else
						fst = -1;
					fst_out << '\t' << fst;
					if (fst > 1)
						cout << "\nFst > 1";
				}
			}
		}
	}
	fst_out.close();


	cout << "\nDone! Input integer to quit\n";
	cin >> end;
	return 0;
}