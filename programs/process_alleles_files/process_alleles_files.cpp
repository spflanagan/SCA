//Author: Sarah P. Flanagan
//Date: 18 September 2015
//Purpose: Take Stacks **.alleles.tsv files and reorient to produce a file with columns:
//IndID, LocusID, Haplotype1, Count1, Haplotype2, Count2

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

int main(int argc, char* argv[])
{
	int end, temp1,temp2;
	int ID, count1, last_loc;
	double temp;
	string hap;
	string alleles_in_name, alleles_out_name, line, whitelist_name;
	ifstream alleles_in, whitelist_file;
	ofstream alleles_out;
	bool interactivemode = false;
	bool whitelist = false;
	string query;
	string tempstring1, tempstring2;
	vector <int> whitelisted_loci;


	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nProcess Alleles Files from Stacks for selection components analysis\n";
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
			cout << "\nProcess Alleles Files from Stacks for selection components analysis\n";
			cout << "-i:\tinput file (with path)\n";
			cout << "-o:\toutput file name (with path)\n";
			cout << "-w:\tOptional whitelist name (with path)\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	alleles_in_name = "default";
	alleles_out_name = "default";

	for (int i = 1; i < argc - 1; i++)
	{
		tempstring1 = argv[i];
		tempstring2 = argv[i + 1];
		if (tempstring1 == "-i")
			alleles_in_name = tempstring2;
		if (tempstring1 == "-o")
			alleles_out_name = tempstring2;
		if (tempstring1 == "-w")
		{
			whitelist = true;
			whitelist_name = tempstring2;
		}

	}

	if (interactivemode)
	{
		cout << "Input Filename:\n";
		cin >> alleles_in_name;
		cout << "\nOutput Filename:\n";
		cin >> alleles_out_name;
		cout << "\nWould you like to use a whitelist?\n";
		cin >> tempstring1;
		if (tempstring1 == "YES" || tempstring1 == "Y" || tempstring1 == "yes" || tempstring1 == "y" || tempstring1 == "Yes")
		{
			whitelist = true;
			cout << "\nWhat is the file name (with path)?\n";
			cin >> whitelist_name;
		}
	}

	if (alleles_in_name == "default")
	{
		cout << "Input Filename:\n";
		cin >> alleles_in_name;
		interactivemode = true;
	}

	if (alleles_out_name == "default")
	{
		cout << "\nOutput Filename:\n";
		cin >> alleles_out_name;
		interactivemode = true;
	}

	cout << "\n\nInput File:\t" << alleles_in_name;
	cout << "\nOutput File:\t" << alleles_out_name;
	if (whitelist)
		cout << "\nUsing whitelist:\t" << whitelist_name;

	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;

		if (query == "n" && query == "N")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> query;
			return 0;
		}
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

	alleles_in.open(alleles_in_name);
	FileTest(alleles_in, alleles_in_name);
	alleles_out.open(alleles_out_name);
	alleles_out << "IndID\tLocID\tHap1\tCount1\tHap2\tCount2";
	last_loc = 0;
	while (!alleles_in.eof())
	{
		while (universal_getline(alleles_in, line))
		{
			stringstream ss;
			ss.str(line);
			ss >> temp1 >> temp2 >> ID >> hap >> temp >> count1;
			if (whitelist)
			{
				bool found = false;
				for (int i = 0; i < whitelisted_loci.size(); i++)
				{
					if (ID == whitelisted_loci[i])
						found = true;
				}
				if (found)
				{
					if (last_loc == ID)
						alleles_out << '\t' << hap << '\t' << count1;
					else
						alleles_out << '\n' << temp2 << '\t' << ID << '\t' << hap << '\t' << count1;
					last_loc = ID;
				}
			}
			else
			{
				if (last_loc == ID)
					alleles_out << '\t' << hap << '\t' << count1;
				else
					alleles_out << '\n' << temp2 << '\t' << ID << '\t' << hap << '\t' << count1;
				last_loc = ID;
			}
		}
	}
	alleles_in.close();
	alleles_out.close();
	cout << "\nOutput file " << alleles_out_name << " closed.\n";
	if (interactivemode)
	{
		cout << "\nFinished! Input integer to exit:\t";
		cin >> end;
		return 0;
	}
	return 0;
}