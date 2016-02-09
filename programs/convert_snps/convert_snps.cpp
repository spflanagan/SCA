//Author: Sarah P. Flanagan
//Date: 25 January 2016
//Purpose: Convert snps file using matches file to a file with the following format:
//SNPID (CatID.Pos), CatalogID, Pos, Allele 1, Allele 1 Count, Allele 2, Allele 2 Count


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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

class locus
{
public:
	int snp_id, cat_id;
	int count;
	string id;
	vector<string> alleles;
	vector<int> allele_count;

	locus()
	{
		count = int();
		cat_id = int();
		snp_id = int();
		id = string();
		alleles = vector < string >();
		allele_count = vector<int>();
	}
};

int main()
{
	int end, temp1, temp2, temp3, temp4, temp5, pos, line_count;
	int ID, count1, loc_index, loc, allele_index;
	size_t t, tt;
	double temp;
	string hap, stemp1, allele, stemp2,stemp3;
	string matches_in_name, snps_in_name, matches_out_name, line;
	ifstream matches_in, snps_in;
	ofstream matches_out;
	vector<locus> individual;


	bool interactivemode = false;
	string query;
	string tempstring1, tempstring2;

	matches_in_name = "E://ubuntushare//SCA//sample_FEM053_align.matches.tsv";
	snps_in_name = "E://ubuntushare//SCA//sample_FEM053_align.snps.tsv";
	matches_out_name = "sample_FEM053.bi.txt";

	snps_in.open(snps_in_name);
	FileTest(snps_in, snps_in_name);
	line_count = 0;
	while (getline(snps_in, line))
	{
		if (!snps_in.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> temp1 >> temp2 >> loc >> pos >> stemp1 >> temp >> hap >> stemp1 >> stemp2 >>stemp3;
			stringstream ssid;
			ssid << loc << "." << pos;
			if (individual.size() > 0)
			{
				loc_index = -5;
				for (t = 0; t < individual.size(); t++)
				{
					
					if (individual[t].id == ssid.str())
						loc_index = t;
				}
				if (loc_index >= 0)
				{
					individual[loc_index].count++;
					allele_index = -5;
					for (tt = 0; tt < individual[loc_index].alleles.size(); tt++)
					{
						if (individual[loc_index].alleles[tt] == hap)
							allele_index = tt;
					}
					if (allele_index < 0)
					{
						individual[loc_index].alleles.push_back(hap);
					}
				}
				else
				{
					individual.push_back(locus());
					individual.back().snp_id = loc;
					individual.back().count = 1;
					individual.back().alleles.push_back(hap);
					individual.back().id = ssid.str();
				}
			}
			else
			{//it's the first one
				individual.push_back(locus());
				individual.back().snp_id = loc;
				individual.back().count = 1;
				individual.back().alleles.push_back(hap);
				individual.back().id = ssid.str();
			}
			line_count++;
			if (line_count % 1000 == 0)
				cout << "\nSuccessfully read in " << line_count << " SNP entries from " << snps_in_name << ".";
		}
	}
	snps_in.close();

	matches_in.open(matches_in_name);
	FileTest(matches_in, matches_in_name);
	line_count = 0;
	while (getline(matches_in, line))
	{
		if (!matches_in.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> temp1 >> temp2 >> loc >> temp3 >> pos >> hap >> count1 >> temp5;
			
			loc_index = -5;
			for (t = 0; t < individual.size(); t++)
			{
				if (individual[t].snp_id == pos)
					individual[t].cat_id = loc;
			}
			if (loc_index >= 0)
			{
				individual[loc_index].count++;
				allele_index = -5;
				for (tt = 0; tt < individual[loc_index].alleles.size(); tt++)
				{
					if (individual[loc_index].alleles[tt] == hap)
						allele_index = tt;
				}
				if (allele_index < 0)
				{
					individual[loc_index].allele_count.push_back(count1);
				}
				else
					individual[loc_index].allele_count[allele_index] = count1;
			}
			line_count++;
			if (line_count % 1000 == 0)
				cout << "\nSuccessfully read in " << line_count << " entries from " << matches_in_name << ".";
		}
	}
	matches_in.close();

	cout << "\nWriting haplotypes to file.\n";
	matches_out.open(matches_out_name);
	matches_out << "SNPID(Cat.Pos)\tAllele1\tAllele1Count\tAllele2\tAllele2Count";
	for (t = 0; t < individual.size(); t++)
	{
		if (individual[t].count == 2)
		{
			matches_out << '\n' << individual[t].id << '\t' << individual[t].alleles[0] << '\t' << individual[t].allele_count[0]
				<< '\t' << individual[t].alleles[1] << '\t' << individual[t].allele_count[1];
		}
		if (individual[t].count == 1)
		{
			matches_out << '\n' << individual[t].id << '\t' << individual[t].alleles[0] << '\t' << individual[t].allele_count[0]
				<< '\t' << individual[t].alleles[0] << '\t' << 0;
		}
	}
	matches_out.close();
	cout << "\nOutput file " << matches_out_name << " closed.\n";
	if (interactivemode)
	{
		cout << "\nFinished! Input integer to exit:\t";
		cin >> end;
		return 0;
	}
	return 0;

}