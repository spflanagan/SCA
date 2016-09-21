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

class fragment
{
public:
	string name;
	int start, end, enz5, enz3;

	fragment()
	{
		name = string();
		start = end = enz5 = enz3 = int();
	}
};


int main()
{
	int i, ii, frag_count,count, start, end,last_enz, this_enz;
	double allelic_dropout_rate, pcr_duplicate_rate, shearing_bias_rate;
	string genome_name, sdigest_name,ddigest_name, line, overhang;
	restriction_enzyme enz1, enz2;
	ifstream genome_file;
	ofstream ddigest_file,sdigest_file;
	vector<fasta_record> genome;
	vector<fragment> digest;
	string sequence, seq_name;

	genome_name = "../../SSC_integrated.fa";
	ddigest_name = "SSC_ddigested.fasta";
	sdigest_name = "SSC_sdigested.fasta";
	enz1.rec_seq = "CTGCAG"; //PstI
	enz1.overhang = "G";
	enz2.rec_seq = "GATC"; //MboI
	enz2.overhang = "";

	//read in the fasta file
	genome_file.open(genome_name);
	FileTest(genome_file, genome_name);
	count = 0;
	sdigest_file.open(sdigest_name);
	sdigest_file << "FragmentID\tFragmentStart\tFragmentEnd\tFragmentLength\t5'Enzyme\t3'Enzyme";
	ddigest_file.open(ddigest_name);
	ddigest_file << "FragmentID\tFragmentStart\tFragmentEnd\tFragmentLength\t5'Enzyme\t3'Enzyme";
	while (!genome_file.eof())
	{		
		universal_getline(genome_file, line);
		if (line.substr(0, 1) == ">" || genome_file.eof())
		{
			if(count > 0)
			{
				//it's the name of the next sequence
				//process previous sequence first
				cout << "\nDigesting " << seq_name << ", which has " << sequence.length() << " characters.";
				start = 0;
				frag_count = 0;
				this_enz = last_enz = 0;
				while (start < sequence.length())
				{
					//double digest
					int first_enz1 = sequence.find(enz1.rec_seq, start);
					int first_enz2 = sequence.find(enz2.rec_seq, start);
					if (first_enz1 < first_enz2)
					{
						end = first_enz1;
						this_enz = 1;
					}
					else
					{
						end = first_enz2;
						this_enz = 2;
					}
					if (end == -1)
					{//stop the loop!
						start = end = sequence.length();
					}
					else
					{
						if (this_enz != last_enz && end - start >= 250 && end - start <= 700)//only keep those with different enzymes
						{																	//and within the desired size range
							stringstream new_name;
							new_name << seq_name << "_frag" << frag_count;
							//record the info in the set of fragments
							digest.push_back(fragment());
							digest.back().name = new_name.str();
							digest.back().start = start;
							digest.back().end = end;
							digest.back().enz5 = last_enz;
							digest.back().enz3 = this_enz;
							//and output it to file
							ddigest_file << '\n' << new_name.str() << '\t' << start << '\t' << end << '\t' << end - start << '\t' << last_enz << '\t' << this_enz;
							//set up for the next time around
							last_enz = this_enz;
							frag_count++;
							if (frag_count % 100000 == 0)
								cout << "\n\tProcessing fragment starting at pos " << start;
						}
						start = end + 1;
					}
						
				}
				cout << "\n\tFound " << frag_count << " fragments on " << seq_name;
			}
				
			//move onto the next one
			count++;
			if (!genome_file.eof())
			{
				seq_name = line.substr(1, line.size());
				sequence.resize(0);
			}

		}
		else
		{//it's the sequence
			if (!genome_file.eof())
			{
				if (line.substr(0, 1) != "\n")
					sequence.append(line);
			}
		}
	}
	genome_file.close();
	ddigest_file.close();
	cout << "\n\nSuccessfully read in " << count << " fasta records.\n";

	//Sample

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}