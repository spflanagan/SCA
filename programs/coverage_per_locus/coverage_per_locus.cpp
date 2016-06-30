//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date: 16 June 2016
//Purpose: Calculate the mean and variance in coverage per locus using matches files from Stacks


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
	int cat_id;
	int count;
	vector<string> alleles;
	vector<int> allele_count;

	locus()
	{
		count = int();
		cat_id = int();
		alleles = vector < string >();
		allele_count = vector<int>();
	}
};

int main()
{
	int end, temp1, temp2, temp3, temp4, temp5;
	int ID, count1, loc_index, loc, allele_index;
	size_t t, tt;
	double temp;
	string hap;
	string matches_in_name, matches_out_name, line, whitelist_name;
	ifstream matches_in, whitelist_file;
	ofstream matches_out;
	vector<locus> individual;
	vector <int> whitelisted_loci;
	bool whitelist = false;
	bool found = false;

	matches_in.open(matches_in_name);
	FileTest(matches_in, matches_in_name);

	while (getline(matches_in, line))
	{
		if (!matches_in.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> temp1 >> temp2 >> loc >> temp3 >> temp4 >> hap >> count1 >> temp5;
			if (whitelist)
			{
				found = false;
				for (int i = 0; i < whitelisted_loci.size(); i++)
				{
					if (whitelisted_loci[i] == loc)
						found = true;
				}
			}
			if (!whitelist || found)
			{
				if (individual.size() > 0)
				{
					loc_index = -5;
					for (t = 0; t < individual.size(); t++)
					{
						if (individual[t].cat_id == loc)
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
							individual[loc_index].allele_count.push_back(count1);
						}
						else
							individual[loc_index].allele_count[allele_index] = count1;
					}
					else
					{
						individual.push_back(locus());
						individual.back().cat_id = loc;
						individual.back().count = 1;
						individual.back().alleles.push_back(hap);
						individual.back().allele_count.push_back(count1);
					}
				}//individual
				else
				{//it's the first one
					individual.push_back(locus());
					individual.back().cat_id = loc;
					individual.back().count = 1;
					individual.back().alleles.push_back(hap);
					individual.back().allele_count.push_back(count1);
				}
			}//whitelist/found
		}
	}
	matches_in.close();

}