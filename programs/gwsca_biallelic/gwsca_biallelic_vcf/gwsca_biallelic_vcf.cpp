//Author: Sarah P. Flanagan
//Date: 10 February 2016
//Purpose: Take reformatted stacks output from many individuals, match it with
//information about their life history stage, and generate Fst values between many 
//different comparisons
//e.g., male-female, adult-offspring, inferred mothers-collected females
//This uses biallelic SNP data gleaned from stacks populations --vcf
//It takes a vcf that includes inferred maternal alleles from infer_mat_vcf
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
#include <algorithm>

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
	vector<int> pop_indices;

	individual()
	{
		filename = ind_id = sex = age = status = string();
		pop_indices = vector<int>();
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
{//need to track a single locus across multiple groups
public:
	string group_name;
	double freq1,freq2;
	double hs, ho, AA,Aa,aa;
	int major_allele_index;
	int two_n, num_inds;
	string allele1, allele2;

	locus_statistics()
	{
		freq1 = freq2 = double();//it's biallelic
		hs = double();
		ho = AA = Aa = aa = double();
		group_name = string();
		two_n = int();
		major_allele_index = int();
		allele1 = allele2 = string();
		num_inds = int();
	}
};

int main()
{
	int i, count, index, pos, line_count, removed_count;
	size_t t, tt, ttt;
	string line, filename, ind_id, snp_id, sex, age, status, allele1, allele2, stemp, quality, chrom;
	string ind_info_name, fst_out_name, summary_out_name, alleles_out_name, vcf_name, debug_out_name;
	ifstream ind_info, vcf;
	ofstream fst_out, summary_out, debug_out;
	vector<int> indices;
	vector<individual> inds;
	vector<locus> reference;
	vector<locus_statistics> pop_stats;

	ind_info_name = "E://ubuntushare//SCA//results//biallelic//ind_info_vcf.txt";//"../../../results/biallelic/id.names.txt";//
	vcf_name = "E://ubuntushare//SCA//results//biallelic//biallelic.gt.vcf";
	summary_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_summary.txt";
	fst_out_name = "E://ubuntushare//SCA//results//biallelic//gwsca_fsts.txt";
	debug_out_name = "../../../results/biallelic/group_assingments.txt";
	
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
				pop_stats[tt].num_inds++;
			}
			if (inds[t].age == pop_stats[tt].group_name)
			{
				age_found = true;
				pop_stats[tt].num_inds++;
			}
			if (inds[t].status == pop_stats[tt].group_name)
			{
				status_found = true;
				pop_stats[tt].num_inds++;
			}
		}
		if (age_found == false)
		{
			if (inds[t].age != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].age;
				pop_stats.back().num_inds = 1;
			}
		}
		if (sex_found == false)
		{
			if (inds[t].sex != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].sex;
				pop_stats.back().num_inds = 1;
			}
		}
		if (status_found == false)
		{
			if (inds[t].status != "0")
			{
				pop_stats.push_back(locus_statistics());
				pop_stats.back().group_name = inds[t].status;
				pop_stats.back().num_inds = 1;
			}
		}
	}
	for (t = 0; t < inds.size(); t++)
	{

		for (tt = 0; tt < pop_stats.size(); tt++)
		{
			if (inds[t].sex == pop_stats[tt].group_name)
				inds[t].pop_indices.push_back(tt);
			if (inds[t].age == pop_stats[tt].group_name)
				inds[t].pop_indices.push_back(tt);
			if (inds[t].status == pop_stats[tt].group_name)
				inds[t].pop_indices.push_back(tt);
		}
	}

	for (i = 0; i < pop_stats.size(); i++)
	{
		pop_stats[i].freq1 = pop_stats[i].freq2 = 0;
		pop_stats[i].allele1 = "0";
		pop_stats[i].allele2 = "1";
		pop_stats[i].two_n = 0;
		pop_stats[i].ho = pop_stats[i].hs = 0;
		pop_stats[i].AA = pop_stats[i].Aa = pop_stats[i].aa = 0;
	}
	cout << "\nEstablished " << pop_stats.size() << " groups:";
	for (t = 0; t < pop_stats.size(); t++)
	{
		cout << '\n' << pop_stats[t].group_name;
	}
	cout << "\n\n";

	fst_out.open(fst_out_name);
	fst_out << "Chrom\tPos";
	for (t = 0; t < pop_stats.size() - 1; t++)
	{
		for (tt = t + 1; tt < pop_stats.size(); tt++)
		{
			fst_out << '\t' << pop_stats[t].group_name << "-" << pop_stats[tt].group_name;
		}
	}

	summary_out.open(summary_out_name);
	summary_out << "Pop\tChrom\tPos\tN\tHs\tHo\tAllele1Freq\tAllele2Freq\tAA\tAa\taa";
	debug_out.open(debug_out_name);
	debug_out << "IndCount\tIndID\tPopIndex\tPopName";
	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	removed_count = line_count = 0;
	while (universal_getline(vcf, line))
	{
		if (!vcf.eof())
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
						bool found = false;
						for (t = 0; t < inds.size(); t++)
						{
							if (stemp == inds[t].filename)
							{
								indices.push_back(t);//then I can just use ind[indices[count]] to access the individuals if necessary.
								found = true;
							}
						}
						if (found == false)
							cout << "WARNING! Individual " << stemp << " not found in " << ind_info_name << " file.\n";
					}
					count++;
				}//while ss >>stemp
			}
			else//then it's a new locus
			{

				stringstream ss;
				ss.str(line);
				ss >> chrom >> pos;
				
				vector<string> alleles;
				count = 0;
				while (ss >> stemp)
				{
					alleles.push_back(stemp);
					if (stemp != "./.")
					{
						find_and_replace(stemp, "/", "\n");
						stringstream ssa;
						ssa.str(stemp);
						getline(ssa, stemp);
						allele1 = stemp;
						getline(ssa, stemp);
						allele2 = stemp;

						//this individual is inds[indices[count]], and belongs to pop_indices.size() groups
						//for each of those groups, add the counts
						for (ttt = 0; ttt < inds[indices[count]].pop_indices.size(); ttt++)
						{
							debug_out << '\n' << count << '\t' << inds[indices[count]].ind_id << '\t' << inds[indices[count]].pop_indices[ttt] << '\t' << pop_stats[inds[indices[count]].pop_indices[ttt]].group_name;
							pop_stats[inds[indices[count]].pop_indices[ttt]].two_n++;
							/*if (pop_stats[inds[indices[count]].pop_indices[ttt]].group_name == "NONPREG")
								cout << "pause for nonpreg\n";*/
							//popstats allele 1
							if (allele1 == "0")
								pop_stats[inds[indices[count]].pop_indices[ttt]].freq1++;
							if (allele2 == "0")
								pop_stats[inds[indices[count]].pop_indices[ttt]].freq1++;
							//popstats allele2
							if (allele1 == "1")
								pop_stats[inds[indices[count]].pop_indices[ttt]].freq2++;
							if (allele2 == "1")
								pop_stats[inds[indices[count]].pop_indices[ttt]].freq2++;

							if (allele1 == allele2 && allele1 != ".")
								pop_stats[inds[indices[count]].pop_indices[ttt]].ho++;
							if (allele1 == "0" && allele2 == "0")
								pop_stats[inds[indices[count]].pop_indices[ttt]].AA++;
							if (allele1 == "0" && allele2 == "1")
								pop_stats[inds[indices[count]].pop_indices[ttt]].Aa++;
							if (allele1 == "1" && allele2 == "0")
								pop_stats[inds[indices[count]].pop_indices[ttt]].Aa++;
							if (allele1 == "1" && allele2 == "1")
								pop_stats[inds[indices[count]].pop_indices[ttt]].aa++;
						}//end popstats
					}//end if
					count++;
				}//end of reading in individuals
				debug_out.close();
				for (i = 0; i < pop_stats.size(); i++)
					pop_stats[i].two_n = pop_stats[i].two_n * 2;
				//now calculate frequencies
				for (t = 0; t < pop_stats.size(); t++)
				{
					//cout << "\n\n" << pop_stats[t].group_name;
					if (pop_stats[t].two_n > 0)//sanity check: make sure the population has members
					{
						if (pop_stats[t].freq1 != 0 && pop_stats[t].freq2 != 0)
						{
							pop_stats[t].ho = 1 - (pop_stats[t].ho / (pop_stats[t].two_n / 2));
							pop_stats[t].AA = pop_stats[t].AA / (pop_stats[t].two_n / 2);
							pop_stats[t].Aa = pop_stats[t].Aa / (pop_stats[t].two_n / 2);
							pop_stats[t].aa = pop_stats[t].aa / (pop_stats[t].two_n / 2);
							if (pop_stats[t].freq1 / pop_stats[t].two_n > 1)
								cout << "pause for large freq.\t";
							pop_stats[t].freq1 = pop_stats[t].freq1 / pop_stats[t].two_n;
							pop_stats[t].freq2 = pop_stats[t].freq2 / pop_stats[t].two_n;
							pop_stats[t].hs = 1 - ((pop_stats[t].freq1 * pop_stats[t].freq1) + (pop_stats[t].freq2*pop_stats[t].freq2));
							summary_out << '\n' << pop_stats[t].group_name << '\t' << chrom << '\t' << pos << '\t' << pop_stats[t].two_n << '\t' << pop_stats[t].hs << '\t'
								<< pop_stats[t].ho << '\t' << pop_stats[t].freq1 << '\t' << pop_stats[t].freq2 << '\t' << pop_stats[t].AA << '\t' << pop_stats[t].Aa << '\t' << pop_stats[t].aa;
						}
						else
							summary_out << '\n' << pop_stats[t].group_name << '\t' << chrom << '\t' << pos << '\t' << pop_stats[t].two_n << '\t' << '\t' << '\t' << pop_stats[t].allele1 << '\t' << pop_stats[t].allele2
							<< '\t' << '\t' << '\t';
					}
					else
						summary_out << '\n' << pop_stats[t].group_name << '\t' << chrom << '\t' << pos << '\t' << pop_stats[t].two_n << '\t' << '\t' << '\t' << pop_stats[t].allele1 << '\t' << pop_stats[t].allele2
						<< '\t' << '\t' << '\t';
				}

				vector<double> fsts;
				vector<double> ht;
				for (t = 0; t < pop_stats.size() -1; t++)
				{
					for (tt = t + 1; tt < pop_stats.size(); tt++)
					{
						if (t != tt)
						{
							fsts.push_back(0);
							ht.push_back(1);
						}
					}
				}

				//calculate ht for every comparison.
				index = 0;
				for (t = 0; t < pop_stats.size() - 1; t++)
				{
					for (tt = t + 1; tt < pop_stats.size(); tt++)
					{
						if (t != tt)
						{
							if (pop_stats[t].freq1 != 0 && pop_stats[t].freq2 != 0 && pop_stats[tt].freq1 != 0 && pop_stats[tt].freq2 != 0)
							{
								double avgp, avgq;
								avgp = avgq = 0;
								avgp = (pop_stats[t].freq1 + pop_stats[tt].freq1) / 2;
								avgq = (pop_stats[t].freq2 + pop_stats[tt].freq2) / 2;
								ht[index] = ht[index] - ((avgp * avgp) + (avgq*avgq));
							}
							else
								ht[index] = -1;
							index++;
						}
					}
				}
				fst_out << '\n' << chrom << '\t' << pos;
				index = 0;
				for (t = 0; t < pop_stats.size() - 1; t++)
				{
					for (tt = t + 1; tt < pop_stats.size(); tt++)
					{
						if (t != tt)
						{
							double hs = (pop_stats[t].hs + pop_stats[tt].hs) / 2;
							if (ht[index] > 0)
								fsts[index] = 1 - (pop_stats[t].hs + pop_stats[tt].hs) / (2*ht[index]);
							else
								fsts[index] = -1.0;
							if (fsts[index]>1)
								cout << "fst is greater than 1\n";
							fst_out << '\t' << fsts[index];
							index++;
						}
					}
				}
				//need to clear out the popstats counters
				for (t = 0; t < pop_stats.size(); t++)
				{
					pop_stats[t].freq1 = pop_stats[t].freq2 = 0;
					pop_stats[t].allele1 = "0";
					pop_stats[t].allele2 = "1";
					pop_stats[t].two_n = 0;
					pop_stats[t].ho = pop_stats[t].hs = 0;
					pop_stats[t].AA = pop_stats[t].Aa = pop_stats[t].aa = 0;
				}
				
				line_count++;
				if (line_count % 1000 == 0)
					cout << "\nCompleted " << line_count << " SNP entries.";
			}//if #		
		}//eof
	}//getline
	vcf.close();
	summary_out.close();
	fst_out.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}