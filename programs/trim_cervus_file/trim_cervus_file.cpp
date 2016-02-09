//Author: Sarah P. Flanagan
//Date: 31 January 2016
//Purpose: Take a Cervus genotypes file and split it into multiple files with a certain number of loci in each.

#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

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
	int i, ii, file_count, locus_count, max_number_loci, count, num_files, index;
	string gen_in_name, gen_out_base, line, temp, id, allele1, allele2;
	ifstream gen_in;
	vector<string> gen_out_names;
	vector<string> locus_names;

	gen_in_name = "../../results/parentage/genotypes.txt";
	gen_out_base = "../../results/parentage/genotypes";
	max_number_loci = 3000;

	gen_in.open(gen_in_name);
	FileTest(gen_in, gen_in_name);
	count = 0;
	while(universal_getline(gen_in, line))
	{
		if (!gen_in.eof())
		{
			if (count == 0)
			{
				stringstream ss;
				ss.str(line);
				ss >> id; //first column is individual ID
				while (ss >> temp)
					locus_names.push_back(temp);
				//now make vectors for output
				num_files = ceil(double(locus_names.size()) / (double(max_number_loci) * 2.0));//there are two alleles per locus
				cout << "There are " << locus_names.size() / 2 << " loci in the dataset, which will be written to " << num_files << " files, each with " << max_number_loci << ".\n";
				index = 0;
				for (i = 0; i < num_files; i++)
				{
					stringstream name;
					name << gen_out_base << i << ".txt";
					gen_out_names.push_back(name.str());
					ofstream out;
					out.open(gen_out_names[i]);
					out << "ID";
					for (ii = 0; ii < max_number_loci * 2; ii++)
					{
						out << '\t' << locus_names[index + ii];
					}
					index = index + ii;
					out.close();
					count++;
				}
			}
			else
			{
				//we write each individual to file.
				stringstream ss;
				ss.str(line);
				ss >> id;
				for (i = 0; i < num_files; i++)
				{
					ofstream out;
					out.open(gen_out_names[i], ios::app);
					out << '\n' << id;
					for (ii = 0; ii < max_number_loci; ii++)
					{
						ss >> allele1 >> allele2;
						out << '\t' << allele1 << '\t' << allele2;
					}
					out.close();
					count++;
				}
			}//end else
		}//end if eof
	}//end while getline
	gen_in.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}