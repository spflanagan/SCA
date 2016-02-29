//Author: Sarah P. Flanagan
//Date: 28 October 2015
//Purpose: Convert matches file to a file with the following format:
//CatalogID, Allele 1, Allele 1 Count, Allele 2, Allele 2 Count


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

int main(int argc, char* argv[])
{
	int end, temp1, temp2, temp3,temp4,temp5;
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
	
	bool interactivemode = false;
	string query;
	string tempstring1, tempstring2;
	
	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nProcess Matches Files from Stacks to get haplotypes per catalog locus\n";
			cout << "-i:\tinput file (with path)\n";
			cout << "-o:\toutput file name (with path)\n";
			cout << "-w:\tOptional whitelist name (with path)\n";
			cout << "no arguments:\tinteractive mode\n";
			cout << "Input integer to quit\n";
			cin >> end;
			return 0;
		}
		if (query == "I" || query == "i")
			interactivemode = true;
	}

	if (argc > 1)
	{
		tempstring1 = argv[1];
		if (tempstring1 == "-h")
		{
			cout << "\nProcess Matches Files from Stacks to get haplotypes per catalog locus\n";
			cout << "-i:\tinput file (with path)\n";
			cout << "-o:\toutput file name (with path)\n";
			cout << "-w:\tOptional whitelist name (with path)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	matches_in_name = "default";
	matches_out_name = "default";

	for (int i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-i")
			matches_in_name = tempstring2;
		if (tempstring1 == "-o")
			matches_out_name = tempstring2;
		if (tempstring1 == "-w")
		{
			whitelist = true;
			whitelist_name = tempstring2;
		}
	}

	if (interactivemode)
	{
		cout << "Input Filename:\n";
		cin >> matches_in_name;
		cout << "\nOutput Filename:\n";
		cin >> matches_out_name;
		cout << "\nWould you like to use a whitelist?\n";
		cin >> tempstring1;
		if (tempstring1 == "YES" || tempstring1 == "Y" || tempstring1 == "yes" || tempstring1 == "y" || tempstring1 == "Yes")
		{
			whitelist = true;
			cout << "\nWhat is the file name (with path)?\n";
			cin >> whitelist_name;
		}
	}

	if (matches_in_name == "default")
	{
		cout << "Input Filename:\n";
		cin >> matches_in_name;
		interactivemode = true;
	}

	if (matches_out_name == "default")
	{
		cout << "\nOutput Filename:\n";
		cin >> matches_out_name;
		interactivemode = true;
	}

	cout << "\n\nInput File:\t" << matches_in_name;
	cout << "\nOutput File:\t" << matches_out_name;
	if (whitelist)
		cout << "\nUsing whitelist:\t" << whitelist_name;

	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query == "n" || query == "N" || query == "NO" || query == "no" || query == "No")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> query;
			return 0;
		}
		else
			cout << "\nProceeding...\n";
	}
	else
		cout << "\nProceeding...\n";


	if (whitelist)
	{
		whitelist_file.open(whitelist_name);
		FileTest(whitelist_file, whitelist_name);
		while (universal_getline(whitelist_file, line))
		{
			if (!whitelist_file.eof())
			{
				whitelisted_loci.push_back(atoi(line.c_str()));
			}
		}
		whitelist_file.close();
	}

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

	cout << "\nWriting haplotypes to file.\n";
	matches_out.open(matches_out_name);
	matches_out << "CatID\tAllele1\tAllele1Count\tAllele2\tAllele2Count";
	for (t = 0; t < individual.size(); t++)
	{
		if (individual[t].count == 2)
		{
			matches_out << '\n' << individual[t].cat_id << '\t' << individual[t].alleles[0] << '\t' << individual[t].allele_count[0]
				<< '\t' << individual[t].alleles[1] << '\t' << individual[t].allele_count[1];
		}
		if (individual[t].count == 1)
		{
			matches_out << '\n' << individual[t].cat_id << '\t' << individual[t].alleles[0] << '\t' << individual[t].allele_count[0]
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
}//end female_list.eof