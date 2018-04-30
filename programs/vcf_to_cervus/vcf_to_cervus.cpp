//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date: 16 June 2016
//Purpose: To take the SNP information from a vcf file and make it compatible with CERVUS.

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

class locus_info
{
public:
	string chrom, locusID;
	int bp;
	vector<char> maternal;
	vector<char> paternal;

	locus_info()
	{
		chrom = locusID = string();
		bp = int();
		maternal = paternal = vector<char>();
	}
};


int main()
{
	int i, ii, count;
	string line, category, all1, all2;
	vector<int> loci_ids;
	char alleles[2];
	vector<string> locusIDs;
	string chrom, locusID, temp;
	vector<string> vcf_ind_ids;
	vector<locus_info> vcf_dat;
	string vcf_name, gen_name;
	ifstream vcf;
	ofstream gen_file;

	vcf_name = "../../results/relatedness/drad_miss3.vcf";//"../../results/stacks/populations_MAF5/batch_1.vcf";
	gen_name = "../../results/parentage_biallelic/snp_genotypes.txt";

	//read it once to get a list of locus IDs
	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	count = 0;
	while (!vcf.eof())
	{
		while (universal_getline(vcf, line))
		{
			if (line.substr(0, 1) != "#" && line != "")
			{
				vcf_dat.push_back(locus_info());
				stringstream ss;
				ss.str(line);
				ss >> vcf_dat[count].chrom >> vcf_dat[count].bp >> vcf_dat[count].locusID >> alleles[0] >> alleles[1] >> temp >> temp >> temp >> temp;
				locusIDs.push_back(vcf_dat[count].locusID);
				while (ss >> temp)
				{
					if (temp.substr(0, 3) != "./.")
					{
						string allele1 = temp.substr(0, 1);
						string allele2 = temp.substr(2, 1);
						char al1 = alleles[atoi(allele1.c_str())];
						char al2 = alleles[atoi(allele2.c_str())];
						vcf_dat[count].maternal.push_back(al1);
						vcf_dat[count].paternal.push_back(al2);
					}
					else
					{
						vcf_dat[count].maternal.push_back('0');
						vcf_dat[count].paternal.push_back('0');
					}
				}
				count++;
			}//if it's a data-filled line
			if (line.substr(0, 2) == "#C")
			{
				stringstream ss;
				ss.str(line);
				ss >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
				while (ss >> temp)
					vcf_ind_ids.push_back(temp);
			}
		}//while vcf line
	}//while vcf
	vcf.close();
	cout << "Found " << locusIDs.size() << " locus records in " << vcf_ind_ids.size() << " individuals.\n";
	
	cout << "Writing data to " << gen_name << ".\n";
	//write genetic info to file
	gen_file.open(gen_name);
	gen_file << "ID";
	for (i = 0; i < locusIDs.size(); i++)
		gen_file << '\t' << locusIDs[i] << "A\t" << locusIDs[i] << "B";
	for (i = 0; i < vcf_ind_ids.size(); i++)
	{
		gen_file << '\n' << vcf_ind_ids[i];
		for (ii = 0; ii < vcf_dat.size(); ii++)
		{
			gen_file << '\t' << vcf_dat[ii].maternal[i] << '\t' << vcf_dat[ii].paternal[i];
		}
	}
	gen_file.close();

	cout << "Done! Input integer to quit.\n";
	cin >> i;
	return 0;

}