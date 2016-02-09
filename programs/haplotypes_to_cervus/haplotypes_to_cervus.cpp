//Author: Sarah P. Flanagan
//Date: 27 January 2016
//Purpose: To take the haplotype information inferred from stacks output and make it compatible with CERVUS.

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

class pairs
{
public:
	string dad, kid;

	pairs()
	{
		dad = string();
		kid = string();
	}
};

class locus_info
{
public:
	int cat_id;
	string allele1, allele2;

	locus_info()
	{
		cat_id = int();
		allele1 = string();
		allele2 = string();
	}
};

int main()
{
	int i, ii, iii, count, locus, cnt1, cnt2, index, cat_id;
	string line, dad, kid, file_name, category, line2, dad_id, ind_id, all1, all2;
	vector<int> loci_ids;
	vector<pairs> dad_kid;
	string locus_file_name, file_list_name, dad_kid_name;
	string off_name, gen_name, par_name;
	ifstream locus_file, file_list, dad_kid_file;
	ofstream gen_file, par_file, off_file;

	locus_file_name = "../../results/haplotypes/loci_summary.txt";
	file_list_name = "../../results/haplotypes/haplotype_files.txt";
	dad_kid_name = "../../results/haplotypes/dad.kid.pairs.txt";

	gen_name = "../../results/parentage/genotypes.txt";
	par_name = "../../results/parentage/candidate_parents.txt";
	off_name = "../../results/parentage/offspring.txt";

	locus_file.open(locus_file_name);
	FileTest(locus_file, locus_file_name);
	count = 0;
	while (universal_getline(locus_file, line))
	{
		if (!locus_file.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> locus;
			loci_ids.push_back(locus);
			count++;
		}
	}
	locus_file.close();
	cout << "\nThere are " << count << " loci in this dataset.\n";

	dad_kid_file.open(dad_kid_name);
	FileTest(dad_kid_file, dad_kid_name);
	count = 0;
	while (universal_getline(dad_kid_file, line))
	{
		if (!dad_kid_file.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> dad >> kid;
			dad_kid.push_back(pairs());
			dad_kid.back().dad = dad;
			dad_kid.back().kid = kid;
			count++;
		}
	}
	dad_kid_file.close();
	cout << "There are " << count << " father-offspring pairs in this dataset.\n";

	gen_file.open(gen_name);
	par_file.open(par_name);
	off_file.open(off_name);
	//write the header lines to file
	gen_file << "ID";
	off_file << "OffspringID\tKnownParentID\t";
	for (i = 0; i < loci_ids.size(); i++)
	{
		gen_file << '\t' << loci_ids[i] << "A" << '\t' << loci_ids[i] << "B";
	}
	//begin processing the files
	file_list.open(file_list_name);
	FileTest(file_list, file_list_name);
	count = 0;
	while (universal_getline(file_list, line))
	{
		if (!file_list.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> file_name >> ind_id >> category;
			
			//read in this file info.
			vector<locus_info> loci; //this will be created new each time for each individual
			ifstream file;
			file.open(file_name);
			FileTest(file, file_name);
			count = 0;
			while (universal_getline(file, line2))
			{
				if (!file.eof())
				{
					if (count > 0)//first row is header
					{
						stringstream ssi;
						ssi.str(line2);
						//read in this haplotype info
						ssi >> cat_id >> all1 >> cnt1 >> all2 >> cnt2;
						loci.push_back(locus_info());
						loci.back().cat_id = cat_id;
						loci.back().allele1 = all1;
						loci.back().allele2 = all2;
					}
					count++;
				}
			}
			file.close();
			//write genetic info to file
			gen_file << '\n' << ind_id;
			//output in order of loci_list
			for (i = 0; i < loci_ids.size(); i++)
			{
				index = -5;
				for (ii = 0; ii < loci.size(); ii++)
				{
					if (loci_ids[i] == loci[ii].cat_id)
						index = ii;
				}
				if (index < 0)
				{
					//if this locus isn't found in this individual, need to output blanks
					gen_file << "\t0\t0";
				}
				else
				{
					gen_file << '\t' << loci[index].allele1 << '\t' << loci[index].allele2;
				}
			}
			if (category == "OFFSPRING" || category == "OFF")
			{
				off_file << "\n" << ind_id;
				//need to find matching known parent ID
				for (i = 0; i < dad_kid.size(); i++)
				{
					if (dad_kid[i].kid == ind_id)
						dad_id = dad_kid[i].kid;
				}
				off_file << '\t' << dad_id;
			}
			if (category == "FEMALE" || category == "FEM")
			{
				if (count == 0)
					par_file << ind_id;
				else
					par_file << '\n' << ind_id;
				count++;
			}
		}
	}
	file_list.close();
	gen_file.close();
	off_file.close();
	par_file.close();

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}