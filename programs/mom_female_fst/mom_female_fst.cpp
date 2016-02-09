//Author: Sarah P. Flanagan
//Date: 1 October 2015
//Purpose: To parse female *.matches.tsv files and MOM*.txt files to compare allele frequencies

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

bool FileTest(ifstream& file, string filename)
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
	return true;
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
	double exp_het;
	int count;
	vector<string> alleles;
	vector<double> freq;

	locus()
	{
		count = int();
		cat_id = int();
		exp_het = double();
		alleles = vector < string >();
		freq = vector<double>();
	}
};

class population
{
public:
	int size;
	vector<locus> loci;

	population()
	{
		size = int();
		loci = vector<locus>();
	}
};

int main()
{
	int end, i, ii, count, loc;
	size_t t, tt;
	vector<population> pops;
	vector<locus> ref;
	string female_list_name, mom_list_name, output_name, summ_out_name;
	string line, filename, allele;
	ifstream female_list, mom_list;
	ofstream output, summ_out;
	int loc_index, allele_index;
	bool found;

	female_list_name = "female.list.txt";
	mom_list_name = "mom.list.txt";
	output_name = "fem.mom.fst.txt";
	summ_out_name = "fem.mom.fst.summary.txt";

	//read in female filenames
	female_list.open(female_list_name);
	FileTest(female_list, female_list_name);
	pops.push_back(population());
	pops[0].size = 0;
	while (getline(female_list, filename))
	{
		if (!female_list.eof())
		{
			ifstream file;
			file.open(filename);
			FileTest(file, filename);
			pops[0].size++;
			while (getline(file, line))
			{
				if (!file.eof())
				{
					stringstream ss;
					ss.str(line);
					ss >> i >> ii >> loc >> i >> ii >> allele >> ii >> i;
					if (pops[0].loci.size() > 0)
					{
						loc_index = -5;
						for (t = 0; t < pops[0].loci.size(); t++)
						{
							if (pops[0].loci[t].cat_id == loc)
								loc_index = t;
						}
						if (loc_index >= 0)
						{
							pops[0].loci[loc_index].count++;
							allele_index = -5;
							for (tt = 0; tt < pops[0].loci[loc_index].alleles.size(); tt++)
							{
								if (pops[0].loci[loc_index].alleles[tt] == allele)
									allele_index = tt;
									
							}
							if (allele_index >= 0)
								pops[0].loci[loc_index].freq[allele_index]++;
							else
							{
								pops[0].loci[loc_index].alleles.push_back(allele);
								pops[0].loci[loc_index].freq.push_back(1);
							}
						}
						else
						{
							pops[0].loci.push_back(locus());
							pops[0].loci.back().cat_id = loc;
							pops[0].loci.back().count = 1;
							pops[0].loci.back().alleles.push_back(allele);
							pops[0].loci.back().freq.push_back(1);
						}
					}
					else
					{//it's the first one
						pops[0].loci.push_back(locus());
						pops[0].loci.back().cat_id = loc;
						pops[0].loci.back().count = 1;
						pops[0].loci.back().alleles.push_back(allele);
						pops[0].loci.back().freq.push_back(1);
					}
				}
			}
			file.close();
		}//end female_list.eof
	}
	female_list.close();
	//read in mom filenames
	mom_list.open(mom_list_name);
	FileTest(mom_list, mom_list_name);
	pops.push_back(population());
	while (getline(mom_list, filename))
	{
		if (!mom_list.eof())
		{
			//open first mom file
			ifstream file;
			file.open(filename);
			FileTest(file, filename);
			pops[1].size++;
			count = 0;
			while (getline(file, line))
			{
				if (!file.eof())
				{
					if (count > 0)
					{
						stringstream ss;
						ss.str(line);
						ss >> loc >> allele;
						if (pops[1].loci.size() > 0)
						{
							loc_index = -5;
							for (t = 0; t < pops[1].loci.size(); t++)
							{
								if (pops[1].loci[t].cat_id == loc)
									loc_index = t;
								
							}
							if (loc_index >= 0)
							{
								pops[1].loci[loc_index].count++;
								allele_index = -5;
								for (tt = 0; tt < pops[1].loci[loc_index].alleles.size(); tt++)
								{
									if (pops[1].loci[loc_index].alleles[tt] == allele)
										allele_index = tt;
								}
								if (allele_index >= 0)
									pops[1].loci[loc_index].freq[allele_index]++;
								else
								{
									pops[1].loci[loc_index].alleles.push_back(allele);
									pops[1].loci[loc_index].freq.push_back(1);
								}
							}
							else
							{
								pops[1].loci.push_back(locus());
								pops[1].loci.back().cat_id = loc;
								pops[1].loci.back().count = 1;
								pops[1].loci.back().alleles.push_back(allele);
								pops[1].loci.back().freq.push_back(1);
							}
						}
						else
						{
							pops[1].loci.push_back(locus());
							pops[1].loci.back().cat_id = loc;
							pops[1].loci.back().count = 1;
							pops[1].loci.back().alleles.push_back(allele);
							pops[1].loci.back().freq.push_back(1);
						}
					}
					count++;
				}
			}
			file.close();
		}//end mom list
	}
	mom_list.close();

	//create the reference
	ref.push_back(locus());
	ref.back().exp_het = 1;
	for (i = 0; i < 2; i++)
	{
		for (t = 0; t < pops[i].loci.size(); t++)
		{
			//see if it's in the reference.
			found = false;
			for (ii = 0; ii < ref.size(); ii++)
			{
				if (ref[ii].cat_id == pops[i].loci[t].cat_id)
				{
					found = true;
					ref[ii].count++;
					for (tt = 0; tt < pops[i].loci[t].alleles.size(); tt++)
					{
						bool allele_found = false;
						for (int j = 0; j < ref[ii].alleles.size(); j++)
						{
							if (pops[i].loci[t].alleles[tt] == ref[ii].alleles[j])
								allele_found = true;
						}
						if (allele_found == false)
						{
							ref[ii].alleles.push_back(pops[i].loci[t].alleles[tt]);
							ref[ii].freq.push_back(0);
						}
					}
				}
			}
			if (found == false)
			{
				ref.push_back(locus());
				ref.back().count = 1;
				ref.back().exp_het = 1;
				ref.back().cat_id = pops[i].loci[t].cat_id;
				for (tt = 0; tt < pops[i].loci[t].alleles.size(); tt++)
				{
					ref.back().alleles.push_back(pops[i].loci[t].alleles[tt]);
					ref.back().freq.push_back(0);
				}
			}
		}
	}
	
	//calculate the frequencies
	summ_out.open(summ_out_name);
	summ_out << "Pop\tLocus\tAllele\tAlleleCount\tFreq\tCount\tPopSize";
	int ref_index;
	for (i = 0; i < 2; i++)
	{
		for (t = 0; t < pops[i].loci.size(); t++)
		{
			pops[i].loci[t].exp_het = 1;
			for (tt = 0; tt < pops[i].loci[t].freq.size(); tt++)
			{
				summ_out << '\n' << i << '\t' << pops[i].loci[t].cat_id << '\t' << pops[i].loci[t].alleles[tt] << '\t' << pops[i].loci[t].freq[tt];
				if (i == 0)
					pops[i].loci[t].freq[tt] = pops[i].loci[t].freq[tt] / (pops[i].size * 2);
				else
					pops[i].loci[t].freq[tt] = pops[i].loci[t].freq[tt] / (pops[i].size);//because moms are essentially haploid
				pops[i].loci[t].exp_het = pops[i].loci[t].exp_het - (pops[i].loci[t].freq[tt] * pops[i].loci[t].freq[tt]);
				//ref[ref_index].exp_het = ref[ref_index].exp_het - (pops[i].loci[t].freq[tt] * pops[i].loci[t].freq[tt]);
				summ_out << '\t' << pops[i].loci[t].freq[tt] <<'\t' << pops[i].loci[t].count << '\t' << pops[i].size;
			}
		}
	}
	summ_out.close();

	//calculate Ht
	int j, jj;
	for (i = 0; i < 2; i++)
	{
		for (t = 0; t < pops[i].loci.size(); t++)
		{
			for (ii = 0; ii < ref.size(); ii++)
			{
				//see if it's in the reference.
				if (ref[ii].cat_id == pops[i].loci[t].cat_id)
				{
					for (tt = 0; tt < pops[i].loci[t].alleles.size(); tt++)
					{
						for (j = 0; j < ref[ii].alleles.size(); j++)
						{
							if (pops[i].loci[t].alleles[tt] == ref[ii].alleles[j])
								ref[ii].freq[j] = ref[ii].freq[j] + pops[i].loci[t].freq[tt];
						}
					}
				}
			}
		}
	}
	for (ii = 0; ii < ref.size(); ii++)
	{		
		for (j = 0; j < ref[ii].alleles.size(); j++)
		{
			ref[ii].freq[j] = ref[ii].freq[j] / 2;
			ref[ii].exp_het = ref[ii].exp_het - (ref[ii].freq[j] * ref[ii].freq[j]);
		}
	}

	double hs;
	cout << "Outputting Fsts\n";
	output.open(output_name);
	output << "CatID\tFst\tHs\tHt";
	vector<locus> fsts;
	for (t = 0; t < pops[0].loci.size(); t++)
	{
		for (tt = 0; tt < pops[1].loci.size(); tt++)
		{
			if (pops[0].loci[t].cat_id == pops[1].loci[tt].cat_id)
			{
				for (ii = 0; ii < ref.size(); ii++)
				{
					if (ref[ii].cat_id == pops[0].loci[t].cat_id || ref[ii].cat_id == pops[1].loci[tt].cat_id)
						ref_index = ii;
				}
				fsts.push_back(locus());
				fsts.back().cat_id = ref[ref_index].cat_id;
				hs = (pops[0].loci[t].exp_het + pops[1].loci[tt].exp_het)/2;
				fsts.back().exp_het = (ref[ref_index].exp_het - hs) / ref[ref_index].exp_het;
				output << '\n' << fsts.back().cat_id << '\t' << fsts.back().exp_het << '\t' << hs << '\t' << ref[ref_index].exp_het;
			}
		}
	}
	output.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}