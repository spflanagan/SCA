//Author: Sarah P. Flanagan
//Last Updated: 13 September 2016
//Date Started: 13 September 2016
//Purpose: Use restriction enzyme recognition sites to digest a genome file to compare sdRAD and ddRAD
//also impose allelic dropout, PCR bias, shearing bias, etc.

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>


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

class restriction_enzyme
{
public:
	string rec_seq;
	string overhang;

	restriction_enzyme()
	{
		overhang = rec_seq = string();
	}

};


int main()
{
	int i, ii, count, start, end;
	string genome_name, digest_name, line, overhang;
	restriction_enzyme enz1, enz2;
	ifstream genome_file;
	ofstream digest_file;
	vector<fasta_record> genome;
	vector<fasta_record> digest;

	genome_name = "../../SSC_integrated.fa";
	digest_name = "SSC_digested.fasta";
	enz1.rec_seq = "CTGCAG"; //PstI
	enz1.overhang = "G";
	enz2.rec_seq = "GATC"; //MboI
	enz2.overhang = "";

	//read in the fasta file
	genome_file.open(genome_name);
	FileTest(genome_file, genome_name);
	count = 0;
	while (!genome_file.eof())
	{
		if (!genome_file.eof())
		{
			universal_getline(genome_file, line);
			if (line.substr(0, 1) == ">")
			{
				//it's the name
				count++;
				genome.push_back(fasta_record());
				genome.back().seq_id = line.substr(1, line.size());
			}
			else
			{//it's the sequence
				if (line.substr(0, 1) != "\n")
					genome.back().sequence.append(line);
			}
		}
	}
	genome_file.close();
	cout << "\nSuccessfully read in " << genome.size() << " fasta records.\n";

	//search through each chrom/scaffold to find restriction sites
	cout << "\nProceeding with in silico digestion, with restriction sites " << enz1.rec_seq << " and " << enz2.rec_seq << '\n';
	digest_file.open(digest_name);
	for (i = 0; i < genome.size(); i++)
	{
		cout << "\n\tDigesting " << genome[i].seq_id;
		start = 0;
		count = 0;
		while (start != genome[i].sequence.length())
		{
			int first_enz1 = genome[i].sequence.find(enz1.rec_seq, start);
			int first_enz2 = genome[i].sequence.find(enz2.rec_seq, start);
			if (first_enz1 < first_enz2)
			{
				end = first_enz1;
				overhang = enz1.overhang;
			}
			else
			{
				end = first_enz2;
				overhang = enz2.overhang;
			}
			//digest.push_back(fasta_record());
			//digest.back().sequence = genome[i].sequence.substr(start, end) + overhang;
			stringstream new_name;
			new_name << genome[i].seq_id << "_frag" << count;
			//digest.back().seq_id = new_name.str();
			digest_file << "\n>" << new_name.str() << '\n' << end - start + 1;// genome[i].sequence.substr(start, end) + overhang;
			count++;
			start = start + 1;
		}
		cout << "\nFound " << count << " fragments on " << genome[i].seq_id;
	}
	digest_file.close();

	cout << "Done! Input integer to quit.\n";
	cin >> i;
	return 0;
}