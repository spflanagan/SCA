//Author: Sarah P. Flanagan
//Date: 24 September 2015
//Purpose: subset a plink file based on a list of individual ID names

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
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

int main()
{
	int end;
	string in_plink_name, list_file_name, out_plink_name, line;
	ifstream in_plink, list_file;
	ofstream out_plink;
	vector<string> list;
	string fam, ind;

	in_plink_name = "E://ubuntushare//SCA//fathers_offspring//batch_1.plink.ped";
	list_file_name = "E://ubuntushare//SCA//fem_names.txt";
	out_plink_name = "E://ubuntushare//SCA//fathers_offspring//females.ped";

	list_file.open(list_file_name);
	FileTest(list_file, list_file_name);
	while (universal_getline(list_file, line))
	{
		if (!list_file.eof())
		{
			list.push_back(line);
		}
	}
	list_file.close();
	cout << "\nList file successfully read! Found " << list.size() << " IDs.\n";

	in_plink.open(in_plink_name);
	FileTest(in_plink, in_plink_name);
	out_plink.open(out_plink_name);
	out_plink << "#Females from SCA";
	while (universal_getline(in_plink, line))
	{
		if (!in_plink.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> fam >> ind;
			bool found = false;
			for (size_t t = 0; t < list.size(); t++)
			{
				if (ind == list[t])
					found = true;
			}
			if (found)
			{
				out_plink << '\n' << line;
				cout << "\nIndividual " << ind << " matched.";
			}
		}
	}
	in_plink.close();
	out_plink.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}