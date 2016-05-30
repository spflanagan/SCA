//Author: Sarah P. Flanagan
//Date: 11 May 2016
//Last updated: 11 May 2016
//Purpose: Take a vcf output from Stacks' populations module and convert to the gentoype format for SNPstats1 (Hohenlohe et al 2010, http://webpages.uidaho.edu/hohenlohe/software.html)


#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>

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
	string chrom, locusID;
	int bp;
	vector<char> maternal;
	vector<char> paternal;

	locus_info()
	{
		chrom = locusID = string();
		bp = int();
		maternal = paternal = vector<char>();
	}
};

class ind_info
{
public:
	string ID, pop;
	int vcf_location;

	ind_info()
	{
		ID = pop = string();
		vcf_location = int();
	}
};

int main()
{
	int i, ii, iii, count;
	char alleles[2];
	string vcf_name, gcall_name, line, temp, popmap_name, popID_name;
	ifstream vcf, popmap;
	ofstream gcall, popID;
	vector<locus_info> vcf_dat;
	vector<ind_info> individuals;
	vector<string> vcf_ind_ids;
	vector<string> pops;

	vcf_name = "../../results/stacks/batch_1.vcf";
	popmap_name = "../../male_female_popmap.txt";
	gcall_name = "../../results/sexlinked/genotype_calls.txt";
	popID_name = "../../results/sexlinked/populationID.txt";

	
	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	count = 0;
	while (!vcf.eof())
	{
		while (universal_getline(vcf, line))
		{
			if (line.substr(0, 1) != "#" && line != "")
			{
				vcf_dat.push_back(locus_info());
				stringstream ss;
				ss.str(line);
				ss >> vcf_dat[count].chrom >> vcf_dat[count].bp >> vcf_dat[count].locusID >> alleles[0] >> alleles[1] >> temp >> temp >> temp >> temp;
				while (ss >> temp)
				{
					if (temp.substr(0, 3) != "./.")
					{
						string allele1 = temp.substr(0, 1);
						string allele2 = temp.substr(2, 1);
						char al1 = alleles[atoi(allele1.c_str())];
						char al2 = alleles[atoi(allele2.c_str())];
						vcf_dat[count].maternal.push_back(al1);
						vcf_dat[count].paternal.push_back(al2);
					}
					else
					{
						vcf_dat[count].maternal.push_back(0);
						vcf_dat[count].paternal.push_back(0);
					}
				}
				count++;
			}//if it's a data-filled line
			if (line.substr(0, 2) == "#C")
			{
				stringstream ss;
				ss.str(line);
				ss >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
				while (ss >> temp)
					vcf_ind_ids.push_back(temp);
			}
		}//while vcf line
	}//while vcf
	vcf.close();
	cout << "Found " << count << " locus records.\n";
	
	popmap.open(popmap_name);
	FileTest(popmap, popmap_name);
	while (!popmap.eof())
	{
		while (universal_getline(popmap, line))
		{
			if (line != "")
			{
				individuals.push_back(ind_info());
				stringstream ss;
				ss.str(line);
				ss >> individuals.back().ID >> individuals.back().pop;
				for (i = 0; i < vcf_ind_ids.size(); i++)
				{
					if (vcf_ind_ids[i] == individuals.back().ID)
						individuals.back().vcf_location = i;
				}
			}
		}
	}
	popmap.close();
	cout << "\nRead in " << individuals.size() << " individual population assignments.\n";
	bool found;
	vector<int> pop_count;
	for (i = 0; i < individuals.size(); i++)
	{
		found = false;
		for (ii = 0; ii < pops.size(); ii++)
		{
			if (individuals[i].pop == pops[ii])
			{
				found = true;
				pop_count[ii]++;
			}
		}
		if (!found)
		{
			pops.push_back(individuals[i].pop);
			pop_count.push_back(1);
		}
	}
	cout << "The individuals belong to " << pops.size() << " populations.\n";

	popID.open(popID_name);
	popID << pops.size();
	int count1;
	for (i = 0; i < pop_count.size(); i++)
		popID << '\n' << pop_count[i];
	for (i = 0; i < pops.size(); i++)
	{
		count1 = 0;
		for (ii = 0; ii < individuals.size(); ii++)
		{
			if (individuals[ii].pop == pops[i])
			{
				if (count1 == 0)
					popID << '\n' << ii;
				else
					popID << '\t' << ii;
				count1++;
			}
		}
	}
	popID.close();

	gcall.open(gcall_name);
	gcall << count << '\t' << individuals.size();
	for (i = 0; i < vcf_dat.size(); i++)
	{
		gcall << '\n' << vcf_dat[i].locusID << '\t' << vcf_dat[i].chrom << '\t' << vcf_dat[i].bp;
		for (ii = 0; ii < individuals.size(); ii++)
		{
			int nt_id, ind_loc;
			ind_loc = individuals[ii].vcf_location;
			if (vcf_dat[i].maternal[ind_loc] == 'A' &&  vcf_dat[i].paternal[ind_loc] == 'A')
				nt_id = 1;
			if (vcf_dat[i].maternal[ind_loc] == 'A' &&  vcf_dat[i].paternal[ind_loc] == 'C')
				nt_id = 2;
			if (vcf_dat[i].maternal[ind_loc] == 'A' &&  vcf_dat[i].paternal[ind_loc] == 'G')
				nt_id = 3;
			if (vcf_dat[i].maternal[ind_loc] == 'A' &&  vcf_dat[i].paternal[ind_loc] == 'T')
				nt_id = 4;
			if (vcf_dat[i].maternal[ind_loc] == 'C' &&  vcf_dat[i].paternal[ind_loc] == 'C')
				nt_id = 5;
			if (vcf_dat[i].maternal[ind_loc] == 'C' &&  vcf_dat[i].paternal[ind_loc] == 'G')
				nt_id = 6;
			if (vcf_dat[i].maternal[ind_loc] == 'C' &&  vcf_dat[i].paternal[ind_loc] == 'T')
				nt_id = 7;
			if (vcf_dat[i].maternal[ind_loc] == 'G' &&  vcf_dat[i].paternal[ind_loc] == 'G')
				nt_id = 8;
			if (vcf_dat[i].maternal[ind_loc] == 'G' &&  vcf_dat[i].paternal[ind_loc] == 'T')
				nt_id = 9;
			if (vcf_dat[i].maternal[ind_loc] == 'T' &&  vcf_dat[i].paternal[ind_loc] == 'T')
				nt_id = 10;
			if (vcf_dat[i].maternal[ind_loc] == '0' &&  vcf_dat[i].paternal[ind_loc] == '0')
				nt_id = 0;
			gcall << '\t' << nt_id;
		}
	}
	gcall.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}