//Author: Sarah P. Flanagan
//Date: 24 September 2015
//Purpose: take lists of loci and alleles and turn it into a plink and a map file.

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

class locus
{
public:
	int loc;
	int distance;
	int bp;

	locus()
	{
		loc = int();
		distance = int();
		bp = int();
	}
};
class map
{
public:
	vector<vector<locus>> chromosome;
	vector<string> chrom_names;

	map()
	{
		chromosome = vector<vector<locus>>();
		chrom_names = vector<string>();
	}
};

class genotypes
{
public:
	vector<string> allele1;
	vector<string> allele2;

	genotypes()
	{
		allele1 = vector<string>();
		allele2 = vector<string>();
	}
};

class ped
{
public:
	string fam_id;
	string ind_id;
	string pat_id;
	string mat_id;
	int sex;//1=male,2=female
	int phenotype;
	vector<genotypes> genotype;

	ped()
	{
		fam_id = string();
		ind_id = string();
		pat_id = string();
		mat_id = string();
		sex = int();
		phenotype = int();
		genotype = vector<genotypes>();
	}
};

int main()
{
	int end, line_count, fem_count;
	int chrom_index, dist, bp, iloc;
	string line, last_chrom, chrom, loc, line2, allele;
	map group_map;
	string map_file_name, mom_plink_name, file_list_name, path;
	ifstream map_file, mom, file_list;
	ofstream mom_plink;

	map_file_name = "E://ubuntushare//SCA//populations//batch_1.plink.map";
	file_list_name = "E://ubuntushare//SCA//fathers_offspring//mom_list.txt";
	mom_plink_name = "E://ubuntushare//SCA//fathers_offspring//mom.plink.ped";
	path = "E://ubuntushare//SCA//fathers_offspring//";

	map_file.open(map_file_name);
	FileTest(map_file, map_file_name);
	last_chrom = "";
	line_count = 0;
	while (universal_getline(map_file, line))
	{
		if (!map_file.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				
				stringstream ss;
				ss.str(line);
				ss >> chrom >> loc >> dist >> bp;
				if (line_count == 0)
				{
					group_map.chromosome.push_back(vector<locus>());
					group_map.chrom_names.push_back(chrom);
					chrom_index = 0;
					group_map.chromosome[chrom_index].push_back(locus());
					group_map.chromosome[chrom_index].back().loc = atoi(loc.substr(0,loc.find("_")).c_str());
					group_map.chromosome[chrom_index].back().bp = bp;
					group_map.chromosome[chrom_index].back().distance = dist;
					last_chrom = chrom;
					line_count++;
				}
				else
				{
					if (chrom != last_chrom)
					{
						group_map.chromosome.push_back(vector<locus>());
						group_map.chrom_names.push_back(chrom);
						chrom_index++;
						last_chrom = chrom;
					}
					group_map.chromosome[chrom_index].push_back(locus());
					group_map.chromosome[chrom_index].back().loc = atoi(loc.substr(0, loc.find("_")).c_str());
					group_map.chromosome[chrom_index].back().distance = dist;
					group_map.chromosome[chrom_index].back().bp = bp;
					line_count++;
				}
			}
		}
	}
	map_file.close();
	cout << "\nMap file read.\n";

	file_list.open(file_list_name);
	FileTest(file_list, file_list_name);
	mom_plink.open(mom_plink_name);
	fem_count = 0;
	while (universal_getline(file_list, line))
	{
		if (!file_list.eof())
		{
			string mom_file_name = path + line;
			ifstream mom_file;
			mom_file.open(mom_file_name);
			FileTest(mom_file, mom_file_name);
			cout << "\nReading file " << mom_file_name << '\n';
			ped mom_ped;
			mom_ped.ind_id = line.substr(0, 6);
			for (size_t t = 0; t < group_map.chromosome.size(); t++)
			{
				mom_ped.genotype.push_back(genotypes());
				for (size_t tt = 0; tt < group_map.chromosome[t].size(); tt++)
				{
					mom_ped.genotype[t].allele1.push_back("0");
					mom_ped.genotype[t].allele2.push_back("0");
				}
			}
			line_count = 0;
			while (universal_getline(mom_file, line2))
			{
				if (!mom_file.eof())
				{
					if (line_count == 0)
						line_count++;
					else
					{
						stringstream ss;
						ss.str(line2);
						ss >> iloc >> allele;
						//find the locus in the map
						for (size_t t = 0; t < group_map.chromosome.size(); t++)
						{
							for (size_t tt = 0; tt < group_map.chromosome[t].size(); tt++)
							{
								if (group_map.chromosome[t][tt].loc == iloc)
								{
									mom_ped.genotype[t].allele1[tt] = allele;
									mom_ped.genotype[t].allele2[tt] = allele;
								}
							}
						}
						line_count++;
					}
				}
			}
			mom_file.close();
			cout << mom_file_name << " had " << line_count << " alleles.\n";
			//write mom to plink file
			if (fem_count == 0)
				mom_plink << "-9\t" << mom_ped.ind_id << "\t-9\t-9\t2\t2";
			else
				mom_plink << "\n-9\t" << mom_ped.ind_id << "\t-9\t-9\t2\t2";
			for (size_t t = 0; t < group_map.chromosome.size(); t++)
			{
				for (size_t tt = 0; tt < group_map.chromosome[t].size(); tt++)
				{
					mom_plink << '\t' << mom_ped.genotype[t].allele1[tt] << '\t' << mom_ped.genotype[t].allele2[tt];
				}
			}
			fem_count++;
		}
	}
	mom_plink.close();
	cout << "Created a plink file with all " << fem_count << " individuals.\n";

	cout << "\nDone! Input integer to quit.\n";
	cin >> end;
	return 0;
}