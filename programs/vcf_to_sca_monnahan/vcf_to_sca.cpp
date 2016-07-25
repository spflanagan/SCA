//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date: 5 July 2016
//Purpose: Convert a vcf file to be able to run it through John Kelly's scripts for the Monnahan et al. (2015) genome-wide selection components analysis

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

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

class locus_id
{
public:
	string chrom, locusID;
	int bp;

	locus_id()
	{
		chrom = locusID = string();
		bp = int();
	}
};

class locus_info
{
public:
	locus_id id;
	vector<char> maternal;
	vector<char> paternal;
	vector<double> LRR;
	vector<double> LRA;
	vector<double> LAA;

	locus_info()
	{
		id = locus_id();
		maternal = paternal = vector<char>();
		LRR = LRA = LAA = vector<double>();
	}
};

class parent_off
{
public:
	string dad_id,kid_id;

	parent_off()
	{
		dad_id = kid_id = string();
	}
};

int main()
{
	int i, ii, iii, count, het_count;
	char alleles[2];
	bool use_whitelist, keep, stacks;
	string vcf_name, line,dad_kid_name,whitelist_name, temp, t2, tt2, out_name;
	ifstream vcf, dad_kid,whitelist;
	ofstream output;
	vector<locus_info> vcf_dat;
	vector<locus_id> keep_loci;
	vector<string> vcf_ind_ids;
	vector<string> pops;
	vector<int> locus_depth;
	vector<double> prop_het;
	vector<double> het_depth_ratio;
	vector<parent_off> dad_kid_info;
	vector<bool> keep_adults;

	use_whitelist = false;
	stacks = true;
	vcf_name = "../../results/stacks/batch_1.vcf";
	dad_kid_name = "../../results/dad.kid.pairs.fullnames.txt";
	out_name = "../../results/monnahan/pipefish_input_full.LLm.txt";
	whitelist_name = "../../results/monnahan/test_whitelist.txt";

	dad_kid.open(dad_kid_name);
	FileTest(dad_kid, dad_kid_name);
	while (!dad_kid.eof())
	{
		while (universal_getline(dad_kid, line))
		{
			if (line != "")
			{
				stringstream ss;
				ss.str(line);
				ss >> t2 >> tt2;
				dad_kid_info.push_back(parent_off());
				dad_kid_info.back().dad_id = t2;
				dad_kid_info.back().kid_id = tt2;
			}
		}
	}
	if (use_whitelist)
	{
		whitelist.open(whitelist_name);
		FileTest(whitelist, whitelist_name);
		while (!whitelist.eof())
		{
			while (universal_getline(whitelist, line))
			{
				if (line != "")
				{
					stringstream ss;
					ss.str(line);
					ss >> t2 >> tt2 >> temp;
					keep_loci.push_back(locus_id());
					keep_loci.back().chrom = t2;
					keep_loci.back().locusID = tt2;
					keep_loci.back().bp = stoi(temp);
				}
			}
		}
		whitelist.close();
		cout << keep_loci.size() << " loci will be retained.";
	}
	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	count = 0;
	while (!vcf.eof())
	{
		while (universal_getline(vcf, line))
		{
			if (line.substr(0, 1) != "#" && line != "")
			{
				
				stringstream ss;
				ss.str(line);
				ss >> t2 >> tt2 >> temp;
				if (use_whitelist)
				{
					keep = false;
					for (i = 0; i < keep_loci.size(); i++)
					{
						if (t2 == keep_loci[i].chrom && stoi(tt2) == keep_loci[i].bp && temp == keep_loci[i].locusID)
							keep = true;
					}
				}
				else
					keep = true;
				if (keep)
				{
					vcf_dat.push_back(locus_info());
					vcf_dat[count].id.chrom = t2;
					vcf_dat[count].id.bp = stoi(tt2);
					vcf_dat[count].id.locusID = temp;
					ss >> alleles[0] >> alleles[1] >> temp >> temp >> temp >> temp;
					while (ss >> temp)
					{
						if (temp.substr(0, 3) != "./.")
						{
							string allele1 = temp.substr(0, 1);
							string allele2 = temp.substr(2, 1);
							stringstream t;
							t.str(temp);
							double d1, d2, d3; 
							if (stacks)
							{
								d1 = d2 = d3 = 0.000001;
								if (allele1 == "0" && allele2 == "0")
									d1 = 1;
								if (allele1 != allele2)
									d2 = 1;
								if (allele1 == "1" && allele2 == "1")
									d3 = 1;
							}
							else
							{
								//if the vcf actually was useful...stacks' output isn't.
								for (i = 0; i < 4; i++)//fourth one is GL
									getline(t, t2, ':');
								stringstream tt;
								tt.str(t2);
								getline(tt, tt2, ',');
								if (tt2 != ".")
									d1 = stod(tt2);
								else d1 = 0.000001;
								getline(tt, tt2, ',');
								if (tt2 != ".")
									d2 = stod(tt2);
								else d2 = 0.000001;
								getline(tt, tt2, ',');
								if (tt2 != ".")
									d3 = stod(tt2);
								else d3 = 0.000001;
							}
							char al1 = alleles[atoi(allele1.c_str())];
							char al2 = alleles[atoi(allele2.c_str())];
							vcf_dat[count].maternal.push_back(al1);
							vcf_dat[count].paternal.push_back(al2);
							vcf_dat[count].LRR.push_back(d1);
							vcf_dat[count].LRA.push_back(d2);
							vcf_dat[count].LAA.push_back(d3);
						}
						else
						{
							vcf_dat[count].maternal.push_back(0);
							vcf_dat[count].paternal.push_back(0);
							vcf_dat[count].LRR.push_back(0);
							vcf_dat[count].LRA.push_back(0);
							vcf_dat[count].LAA.push_back(0);
						}
					}
					count++;
				}
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
	cout << "Found " << count << " locus records.\n";

	for (i = 0; i < vcf_ind_ids.size(); i++)
	{
		stringstream ss;
		ss.str(vcf_ind_ids[i]);
		getline(ss, line,'_');
		getline(ss, line,'_');
		if (line.substr(0, 3) == "PRM" || line.substr(0, 3) == "FEM")
			keep_adults.push_back(true);
		else
			keep_adults.push_back(false);
	}
	//write to file
	output.open(out_name);
	cout << "Writing output to " << out_name;
	//chromosome,position,name,bases,number of progeny, mom - LRR,mom - LRA,mom - LAA,sire ID for 1st progeny,offspring1 - LRR, offspring1 - LRA,offspring1 - ​LAA,sire ID for 2st progeny,offspring2 - LRR,offspring2 - LRA,offspring2​ - ​LAA
	for (i = 0; i < vcf_dat.size(); i++)
	{
		for (ii = 0; ii < vcf_ind_ids.size(); ii++)
		{
			if (keep_adults[ii] == true)
			{
				if (vcf_ind_ids[ii].substr(7, 3) == "FEM")
				{//if it's female
					stringstream ss;
					ss.str(vcf_ind_ids[ii]);
					getline(ss, line, '_');
					getline(ss, line, '_');
					if (vcf_dat[i].maternal[ii] != '\0')
						output << '\n' << vcf_dat[i].id.chrom << '\t' << vcf_dat[i].id.bp << '\t' << line.substr(0, 3) << "." << line.substr(3, line.length())
							<< '\t' << vcf_dat[i].maternal[ii] << "," << vcf_dat[i].paternal[ii] << '\t'
							<< vcf_dat[i].LRR[ii] << '\t' << vcf_dat[i].LRA[ii] << '\t' << vcf_dat[i].LAA[ii];
				}
				else//it's male 
				{
					vector<int> kid_id;
					bool skip = false;
					for (iii = 0; iii < dad_kid_info.size(); iii++)//check to see if this ind is a dad or a kid.
					{
						if (vcf_ind_ids[ii] == dad_kid_info[iii].dad_id)
						{
							for (int iv = 0; iv < vcf_ind_ids.size(); iv++)
							{
								if (dad_kid_info[iii].kid_id == vcf_ind_ids[iv] && vcf_dat[i].maternal[iv] != '\0')
									kid_id.push_back(iv);
							}
						}
					}
					if (vcf_dat[i].maternal[ii] != '\0')
					{
						stringstream ss;
						ss.str(vcf_ind_ids[ii]);
						getline(ss, line, '_');
						getline(ss, line, '_');
						output << '\n' << vcf_dat[i].id.chrom << '\t' << vcf_dat[i].id.bp << '\t' << line.substr(0, 3) << "." << line.substr(3, line.length())
							<< '\t' << vcf_dat[i].maternal[ii] << "," << vcf_dat[i].paternal[ii] << '\t';
						if (kid_id.size() > 0)
							output << kid_id.size()<< '\t';
						else
							output << "0\t";
						output << vcf_dat[i].LRR[ii] << '\t' << vcf_dat[i].LRA[ii] << '\t' << vcf_dat[i].LAA[ii];
						for (iii = 0; iii < kid_id.size(); iii++)
						{
							if (vcf_dat[i].maternal[kid_id[iii]] != '\0')
								output << "\t1\t" << vcf_dat[i].LRR[kid_id[iii]] << '\t' << vcf_dat[i].LRA[kid_id[iii]] << '\t' << vcf_dat[i].LAA[kid_id[iii]];
						}
					}
					else
					{//the adult wasn't genotyped
						if (kid_id.size() > 0)
						{
							skip = true;
							for (iii = 0; iii < kid_id.size(); iii++)
							{
								if (vcf_dat[i].maternal[kid_id[iii]] != '\0')
									skip = false;
							}
							if (!skip && kid_id.size() > 0)
							{
								stringstream ss;
								ss.str(vcf_ind_ids[ii]);
								getline(ss, line, '_');
								getline(ss, line, '_');
								output << '\n' << vcf_dat[i].id.chrom << '\t' << vcf_dat[i].id.bp << '\t' << line.substr(0, 3) << "." << line.substr(3, line.length())
									<< '\t' << "X,X" << '\t';
								if (kid_id.size() > 0)
									output << kid_id.size() << '\t';
								else
									output << "0\t";
								output << 1 << '\t' << 1 << '\t' << 1;
								for (iii = 0; iii < kid_id.size(); iii++)
								{
									if (vcf_dat[i].maternal[kid_id[iii]] != '\0')
										output << "\t1\t" << vcf_dat[i].LRR[kid_id[iii]] << '\t' << vcf_dat[i].LRA[kid_id[iii]] << '\t' << vcf_dat[i].LAA[kid_id[iii]];
								}
							}
						}
						
					}
				}
			}
		}
	}
	
	output.close();

	cout << "\nDone!";
	cin >> i;
	return 0;
}