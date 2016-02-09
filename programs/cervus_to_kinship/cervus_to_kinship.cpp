//Author: Sarah P. Flanagan
//Date: 3 February 2016
//Purpose: Take a CERVUS input file and convert it to a KINSHIP input file

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

int main()
{
	int i, ii, count;
	string cervus_name, kinship_name, line, temp1, temp2;
	ifstream cervus_file;
	ofstream kinship_file;

	cervus_name = "../../results/parentage/relatedness1000.txt";
	kinship_name = "../../results/parentage/relatedness1000_kin.txt";

	cervus_file.open(cervus_name);
	FileTest(cervus_file, cervus_name);
	kinship_file.open(kinship_name);
	count = 0;
	while (universal_getline(cervus_file, line))
	{
		if (!cervus_file.eof())
		{
			stringstream ss;
			ss.str(line);
			if (count == 0)
			{
				ss >> temp1;
				kinship_file << temp1;
				while (ss >> temp1 >> temp2)
					kinship_file << '\t' << temp1.substr(0, temp1.size() - 1);
			}
			else
			{
				ss >> temp1;
				kinship_file << '\n' <<  temp1;
				while (ss >> temp1 >> temp2)
				{
					if (temp1 != "0")
						kinship_file << '\t' << temp1 << "/" << temp2;
					else
						kinship_file << "\t/";
				}
			}
			count++;
		}
	}
	cervus_file.close();
	kinship_file.close();
	cout << "Wrote " << count << " lines to " << kinship_name;

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}