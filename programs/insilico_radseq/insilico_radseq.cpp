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
#include <time.h>
#include <random>
#include "random_numbers.h"


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

class restriction_site
{
public:
	bool polymorphic;
	double mutation_rate;
	int position;
	string enzyme;
	bool sd_tag, dd_tag;

	restriction_site()
	{
		polymorphic = sd_tag = dd_tag = bool();
		mutation_rate = double();
		position = int();
		enzyme = string();
	}
};

class fragment
{
public:
	string name;
	int start, end, enz5, enz3, pos;
	string chrom;
	bool polymorphic, shared;
	double mutation_rate;

	fragment()
	{
		name = string();
		pos = start = end = enz5 = enz3 = int();
		chrom = string();
		polymorphic =shared= bool();
		mutation_rate = double();
	}
};

class individual
{
public:
	vector<int> allele1, allele2;

	individual()
	{
		allele1 = allele2 = vector<int>();
	}
};

int main()
{
	int i, ii, frag_count,count, start, end,last_enz, this_enz,Ne;
	int s_start, s_end, s_frag_count, s_cutsite, s_length, sd_ind, dd_ind;
	double mutation_shape_param, pcr_duplicate_rate, shearing_bias_rate, mean_sheared_length;
	double mutation_rate;
	string genome_name, sdigest_name,ddigest_name, line, overhang;
	restriction_enzyme enz1, enz2;
	ifstream genome_file;
	ofstream ddigest_file,sdigest_file;
	vector<fasta_record> genome;
	vector<fragment> ddigest, sdigest;
	vector<individual> sd_population, dd_population;
	vector<restriction_site> digestion;
	string sequence, seq_name;
	sgenrand(time(0));

	
	genome_name = "../../SSC_integrated.fa";
	ddigest_name = "../../results/insilico/SSC_ddigested.txt";
	sdigest_name = "../../results/insilico/SSC_sdigested.txt";
	enz1.rec_seq = "CTGCAG"; //PstI
	enz1.overhang = "G";
	enz2.rec_seq = "GATC"; //MboI
	enz2.overhang = "";
	sd_ind = dd_ind = 50;
	pcr_duplicate_rate = 0.05;
	Ne = 10000;

	//read in the fasta file
	genome_file.open(genome_name);
	FileTest(genome_file, genome_name);
	count = 0;
	sdigest_file.open(sdigest_name);
	sdigest_file << "FragmentID\tFragmentStart\tFragmentEnd\tFragmentLength\t5'Enzyme\t3'Enzyme\tChrom\tPstIPos";
	ddigest_file.open(ddigest_name);
	ddigest_file << "FragmentID\tFragmentStart\tFragmentEnd\tFragmentLength\t5'Enzyme\t3'Enzyme\tChrom\tPstIPos";
	while (!genome_file.eof())
	{		
		universal_getline(genome_file, line);
		if (line.substr(0, 1) == ">" || genome_file.eof())
		{
			if(count > 0)
			{
				//it's the name of the next sequence
				//process previous sequence first
				
				cout << "\nsingle-digesting " << seq_name << ", which has " << sequence.length() << " characters.";
				//single digest
				s_start = s_frag_count = 0;
				mean_sheared_length = 0;
				while (s_start < sequence.length())
				{
					s_cutsite = sequence.find(enz1.rec_seq, s_start);
					if (s_cutsite == -1)
					{
						s_start = s_end = s_cutsite = sequence.length();
					}
					else
					{
						s_length = s_cutsite - s_start;
						if ((mean_sheared_length + s_length)/(s_frag_count+1) > 500)//then it gets sheared
						{															//and we only keep sheared fragments
							int new_length = s_length;
							int new_end = 0;
							int new_start = s_start;
							int num_shears = s_length / 500.0;
							int min_dist, max_dist, shear_site;
							min_dist= s_length;
							max_dist = 0;
							string new_seq1, new_seq,new_seq2;
							new_seq = sequence.substr(s_start,s_cutsite-s_start);
							for (i = 0; i < num_shears; i++)
							{
								shear_site = randnum(new_length);
								if (shear_site < min_dist)
									min_dist = shear_site;
								if (shear_site > max_dist)
									max_dist = shear_site;

								////split the sequence 
								//new_seq1 = new_seq.substr(0, new_end);
								//new_seq2 = new_seq.substr(new_end, new_length - new_end+1);
								//if (i == num_shears - 1)//if it's the last one, add both
								//{
								//	//evaluate seq2
								//	if (s_start == (new_start + new_end) || (new_start + new_length) == s_cutsite)
								//	{
								//		if ((new_length - new_end) >= 250 && (new_length - new_end) <= 700)
								//		{
								//			sdigest.push_back(fragment());
								//			stringstream new_name;
								//			new_name << seq_name << "_sd_frag" << s_frag_count;
								//			sdigest.back().name = new_name.str();
								//			sdigest.back().start = new_start + new_end;
								//			sdigest.back().end = new_start + new_length;
								//			s_frag_count++;
								//			mean_sheared_length = mean_sheared_length + new_seq2.length();
								//			if (s_start == sdigest.back().start)
								//			{
								//				sdigest.back().enz5 = 1;//it's the last cutsite
								//				sdigest.back().pos = s_start;
								//				sdigest.back().chrom = seq_name;
								//			}
								//			else
								//				sdigest.back().enz5 = 2; //it's a sheared end
								//			if (sdigest.back().end == s_cutsite)
								//			{
								//				sdigest.back().enz3 = 1;//it's this cutsite
								//				sdigest.back().pos = s_cutsite;
								//				sdigest.back().chrom = seq_name;
								//			}
								//			else
								//				sdigest.back().enz3 = 2;//it's sheared
								//			sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
								//				<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
								//				<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
								//		}
								//	}
								//	//add seq1
								//	if (s_start == new_start || (new_start + new_end) == s_cutsite)
								//	{
								//		if (new_end >= 250 && new_end <= 700)
								//		{
								//			sdigest.push_back(fragment());
								//			stringstream new_name;
								//			new_name << seq_name << "_sd_frag" << s_frag_count;
								//			sdigest.back().name = new_name.str();
								//			sdigest.back().start = new_start;
								//			sdigest.back().end = new_start + new_end;
								//			s_frag_count++;
								//			mean_sheared_length = mean_sheared_length + new_seq1.length();
								//			if (s_start == sdigest.back().start)
								//			{
								//				sdigest.back().enz5 = 1;//it's the last cutsite
								//				sdigest.back().pos = s_start;
								//				sdigest.back().chrom = seq_name;
								//			}
								//			else
								//				sdigest.back().enz5 = 2; //it's a sheared end
								//			if (sdigest.back().end == s_cutsite)
								//			{
								//				sdigest.back().enz3 = 1;//it's this cutsite
								//				sdigest.back().pos = s_cutsite;
								//				sdigest.back().chrom = seq_name;
								//			}
								//			else
								//				sdigest.back().enz3 = 2;//it's sheared
								//			sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
								//				<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
								//				<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
								//		}
								//	}
								//}
								//else
								//{
								//	if (new_seq1.length() > new_seq2.length())
								//	{
								//		if (s_start == (new_start + new_end) || (new_start + new_length) == s_cutsite)
								//		{
								//			new_seq = new_seq1;//new_seq1 will be further sheared
								//			if ((new_length - new_end) >= 250 && (new_length - new_end) <= 700)
								//			{
								//				sdigest.push_back(fragment());
								//				stringstream new_name;
								//				new_name << seq_name << "_sd_frag" << s_frag_count;
								//				sdigest.back().name = new_name.str();
								//				sdigest.back().start = new_start + new_end;
								//				sdigest.back().end = new_start + new_length;
								//				s_frag_count++;
								//				new_start = new_start;//it doesn't change
								//				mean_sheared_length = mean_sheared_length + new_seq2.length();
								//				if (s_start == sdigest.back().start)
								//				{
								//					sdigest.back().enz5 = 1;//it's the last cutsite
								//					sdigest.back().pos = s_start;
								//					sdigest.back().chrom = seq_name;
								//				}
								//				else
								//					sdigest.back().enz5 = 2; //it's a sheared end
								//				if (sdigest.back().end == s_cutsite)
								//				{
								//					sdigest.back().enz3 = 1;//it's this cutsite
								//					sdigest.back().pos = s_cutsite;
								//					sdigest.back().chrom = seq_name;
								//				}
								//				else
								//					sdigest.back().enz3 = 2;//it's sheared
								//				sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
								//					<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
								//					<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
								//			}
								//		}
								//	}
								//	else
								//	{
								//		if (s_start == new_start ||( new_start + new_end) == s_cutsite)
								//		{
								//			new_seq = new_seq2;//new_seq2 will be further sheared
								//			if (new_end >= 250 && new_end <= 700)
								//			{
								//				sdigest.push_back(fragment());
								//				stringstream new_name;
								//				new_name << seq_name << "_sd_frag" << s_frag_count;
								//				sdigest.back().name = new_name.str();
								//				sdigest.back().start = new_start;
								//				sdigest.back().end = new_start + new_end;
								//				s_frag_count++;
								//				new_start = new_start + new_end + 1;
								//				mean_sheared_length = mean_sheared_length + new_seq1.length();
								//				if (s_start == sdigest.back().start)
								//				{
								//					sdigest.back().enz5 = 1;//it's the last cutsite
								//					sdigest.back().pos = s_start;
								//					sdigest.back().chrom = seq_name;
								//				}
								//				else
								//					sdigest.back().enz5 = 2; //it's a sheared end
								//				if (sdigest.back().end == s_cutsite)
								//				{
								//					sdigest.back().enz3 = 1;//it's this cutsite
								//					sdigest.back().pos = s_cutsite;
								//					sdigest.back().chrom = seq_name;
								//				}
								//				else
								//					sdigest.back().enz3 = 2;//it's sheared
								//				sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
								//					<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
								//					<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
								//			}
								//		}
								//	}
								//}
								//new_length = new_seq.length();
							}//num_shears
							if (min_dist >= 250 && min_dist <= 700)//then we keep the first sheared frag
							{
								sdigest.push_back(fragment());
								stringstream new_name;
								new_name << seq_name << "_sd_frag" << s_frag_count;
								sdigest.back().name = new_name.str();
								sdigest.back().start = s_start;
								sdigest.back().end = s_start + min_dist;
								s_frag_count++;
								mean_sheared_length = mean_sheared_length + min_dist;
								sdigest.back().enz5 = 1;//it's the last cutsite
								sdigest.back().pos = s_start;
								sdigest.back().chrom = seq_name;
								sdigest.back().enz3 = 4;//it's sheared
								sdigest.back().polymorphic = false;
								sdigest.back().shared = false;
								sdigest.back().mutation_rate = 0;
								sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
									<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
									<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
							}
							if ((s_cutsite - max_dist) >= 250 && (s_cutsite - max_dist <= 700))//then we keep the last fragment
							{
								sdigest.push_back(fragment());
								stringstream new_name;
								new_name << seq_name << "_sd_frag" << s_frag_count;
								sdigest.back().name = new_name.str();
								sdigest.back().start = s_start + max_dist;
								sdigest.back().end = s_cutsite;
								s_frag_count++;
								mean_sheared_length = mean_sheared_length + min_dist;
								sdigest.back().enz5 = 4;//it's the last cutsite
								sdigest.back().pos = s_start + max_dist;
								sdigest.back().chrom = seq_name;
								sdigest.back().enz3 = 1;//it's sheared
								sdigest.back().polymorphic = false;
								sdigest.back().shared = false;
								sdigest.back().mutation_rate = 0;
								sdigest_file << '\n' << sdigest.back().name << '\t' << sdigest.back().start << '\t' << sdigest.back().end << '\t'
									<< sdigest.back().end - sdigest.back().start << '\t' << sdigest.back().enz5 << '\t' << sdigest.back().enz5
									<< '\t' << sdigest.back().chrom << '\t' << sdigest.back().pos;
							}
						}
						s_start = s_cutsite + 1;
					}
				}//s_start
				cout << "\n\tFound " << s_frag_count << " single-digest fragments on " << seq_name << " with mean fragment size " << mean_sheared_length/s_frag_count;
				//double digest
				cout << "\ndouble-digesting " << seq_name << ", which has " << sequence.length() << " characters.";
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
						if (this_enz != last_enz && (end - start) >= 250 && (end - start) <= 700)//only keep those with different enzymes
						{																	//and within the desired size range
							stringstream new_name;
							new_name << seq_name << "_frag" << frag_count;
							//record the info in the set of fragments
							ddigest.push_back(fragment());
							ddigest.back().name = new_name.str();
							ddigest.back().start = start;
							ddigest.back().end = end;
							ddigest.back().enz5 = last_enz;
							ddigest.back().enz3 = this_enz;
							ddigest.back().shared = false;
							if (this_enz == 1)
								ddigest.back().pos = end;
							else
								ddigest.back().pos = start;
							ddigest.back().chrom = seq_name;
							//and output it to file
							ddigest_file << '\n' << new_name.str() << '\t' << start << '\t' << end << '\t' << end - start << '\t' << last_enz << '\t' << this_enz 
								<< '\t' << ddigest.back().chrom << '\t' << ddigest.back().pos;
							//set up for the next time around
							last_enz = this_enz;
							frag_count++;
							if (frag_count % 100000 == 0)
								cout << "\n\tProcessing fragment starting at pos " << start;
						}
						start = end + 1;
					}
						
				}
				cout << "\n\tFound " << frag_count << " double-digest fragments on " << seq_name;
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
	sdigest_file.close();
	cout << "\n\nSuccessfully read in " << count << " fasta records.\n";

	//set up population counters
	for (i = 0; i < sd_ind; i++)
	{
		sd_population.push_back(individual());		
	}
	for (i = 0; i < dd_ind; i++)
	{
		dd_population.push_back(individual());
	}
	//set up the mutational structure of restriction sites
	for (i = 0; i < ddigest.size(); i++)
	{
		if (genrand() > 0.1)//if not constant, could be null
			ddigest[i].polymorphic = true;
		else
			ddigest[i].polymorphic = false;
		default_random_engine generator;
		uniform_real_distribution<double> distribution(0.000000001, 0.00000001);
		ddigest[i].mutation_rate = 4*Ne*distribution(generator);
		for (ii = 0; ii < sdigest.size(); ii++)
		{
			if (ddigest[i].chrom == sdigest[ii].chrom && ddigest[i].pos == sdigest[ii].pos)
			{
				sdigest[ii].polymorphic = ddigest[i].polymorphic;
				sdigest[ii].mutation_rate = ddigest[i].mutation_rate;
				sdigest[ii].shared = true;
				ddigest[i].shared = true;
			}
		}
	}
	//now fill in the blanks in the sdigest
	for (i = 0; i < sdigest.size(); i++)
	{
		if (sdigest[i].shared == false)
		{
			if (genrand() > 0.1)//if not constant, could be null
				sdigest[i].polymorphic = true;
			else
				sdigest[i].polymorphic = false;
			default_random_engine generator;
			uniform_real_distribution<double> distribution(0.000000001, 0.00000001);
			sdigest[i].mutation_rate = 4 * Ne*distribution(generator);
		}
	}
	//now time to sample (include pcr bias?)

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}