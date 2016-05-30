//Author: Sarah P. Flanagan
//Date: 2 November 2015
//Purpose: Extract part of a fasta sequence from a starting bp to an ending bp

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


class fasta_record
{
public:
	string seq_id;
	string sequence;

	fasta_record()
	{
		seq_id = "";
		sequence = "";
	}
};

class subset_info
{
public:
	string chrom_name;
	int start, end;

	subset_info()
	{
		chrom_name = string();
		start = end = int();
	}
};

int main(int argc, char* argv[])
{
	int end, start, count, i, ii, iii;
	string fasta_name, line, fasta_subset_name, seq_name, input_name;
	ifstream fasta, input;
	ofstream fasta_subset;
	vector<fasta_record> seq;
	vector<subset_info> subsetting;
	string tempstring1, tempstring2, query;
	bool interactivemode;

	if (argc == 1)
	{
		cout << "\n(I)nteractive or (H)elp?\n";
		cin >> query;
		if (query == "H" || query == "h")
		{
			cout << "\nextract_sequence_part\n";
			cout << "Extract from a fasta sequence a specific region in BP number\n";
			cout << "-f:\tfasta file input\n";
			cout << "-i:\tInput file\n";
			cout << "-o:\tOutput file name.\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
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
			cout << "\nextract_sequence_part\n";
			cout << "Extract from a fasta sequence a specific region in BP number\n";
			cout << "-f:\tfasta file input\n";
			cout << "-i:\tInput file\n";
			cout << "-o:\tOutput file name.\n";
			cout << "-h:\tPrints this message\n";
			cout << "no arguments:\tinteractive mode\n";
			return 0;
		}
	}

	for (i = 1; i < argc - 1; i++)
	{
		interactivemode = false;
		tempstring1 = argv[i];
		//tempstring2 = argv[i + 1];
		if (tempstring1 == "-f")
			fasta_name = argv[i + 1];
		if(tempstring1 == "-i")
			input_name = argv[i + 1];
		if (tempstring1 == "-o")
			fasta_subset_name = argv[i + 1];
	}

	if (interactivemode)
	{
		cout << "\nEnter file with a fasta seqence.\n";
		cin >> fasta_name;
		cout << "\nEnter input file name.\n";
		cin >> input_name;
		cout << "\nProvide output file name\n";
		cin >> fasta_subset_name;
	}

	cout << "\nFasta file:\t" << fasta_name;
	cout << "\nInput name:\t" << input_name;
	cout << "\nOutput name:\t" << fasta_subset_name;
	if (interactivemode)
	{
		cout << "\n\nProceed? (y to proceed)\n";
		cin >> query;
		if (query != "y" && query != "Y")
		{
			cout << "\n\nEnter an integer to exit!!\n";
			cin >> i;
			return 0;
		}
	}
	else
		cout << "\n\nProceeding...\n";
	
	input.open(input_name);
	FileTest(input, input_name);
	while (!input.eof())
	{
		while (universal_getline(input, line))
		{
			if (line != "")
			{
				subsetting.push_back(subset_info());
				stringstream ss;
				ss.str(line);
				ss >> subsetting.back().chrom_name >> subsetting.back().start >> subsetting.back().end;
			}
		}
	}
	input.close();
	cout << "\nFound " << subsetting.size() << " regions to extract.\n";

	fasta.open(fasta_name);
	FileTest(fasta, fasta_name);
	count = 0;
	while (!fasta.eof())
	{
		if (!fasta.eof())
		{
			universal_getline(fasta, line);
			if (line.substr(0, 1) == ">")
			{
				//it's the name
				count++;
				seq.push_back(fasta_record());
				seq.back().seq_id = line.substr(1, line.size());
			}
			else
			{//it's the sequence
				if (line.substr(0, 1) != "\n")
					seq.back().sequence.append(line);
			}
		}
	}
	fasta.close();

	
	fasta_subset.open(fasta_subset_name);
	count = 0;
	for (iii = 0; iii < subsetting.size(); iii++)
	{
		for (ii = 0; ii < seq.size(); ii++)
		{
			if (seq[ii].seq_id == subsetting[iii].chrom_name)
			{
				if (count ==0)
					fasta_subset << ">" << seq[ii].seq_id << "_" << subsetting[iii].start << "-" << subsetting[iii].end << '\n';
				else
					fasta_subset << "\n>" << seq[ii].seq_id << "_" << subsetting[iii].start << "-" << subsetting[iii].end << '\n';
				for (i = 0; i < seq[ii].sequence.length(); i++)
				{
					if (i > subsetting[iii].start && i < subsetting[iii].end)
						fasta_subset << seq[ii].sequence[i];
				}
				count++;
			}
		}
	}
	fasta_subset.close();

	cout << "\nDone!";
	if (interactivemode)
	{
		cout << "Enter integer to quit.\n";
		cin >> end;
		return 0;
	}
	else
		return 0;
}