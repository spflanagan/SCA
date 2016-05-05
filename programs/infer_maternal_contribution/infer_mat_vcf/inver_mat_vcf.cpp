//Date: 18 January 2016
//Author: Sarah P. Flanagan (sflanagan@bio.tamu.edu)
//Purpose: Using a vcf file, infer maternal phenotype

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
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

class locus_data
{
public:
	int pos, id; 
	string chrom, ref, alt;
	
	locus_data()
	{
		chrom = string();
		ref = string();
		alt = string();
		pos = int();
		id = int();
	}
};

class locus
{
public:
	int index;
	string allele;

	locus()
	{
		index = int();
		allele = string();
	}
};

class comparisons
{
public:
	int dad_index, kid_index;
	string dad_name, kid_name;
	vector <string> mom_allele; //used to be locus

	comparisons()
	{
		dad_index = int();
		kid_index = int();
		dad_name = string();
		kid_name = string();
		mom_allele = vector<string>();
	}
};



int main()
{
	int i, ii, iii;
	int num_dads, num_off, gt_index, locus_index;
	int line_count, snp_count, count, prog_start, prog_end;
	string vcf_name, output_name, line, stemp;
	string chrom, ID, ref, alt, filter, info, format, temp, temp2, last_chr;
	int pos, index, num_removed = 0;
	string quality, dad_kid_name, allele1, allele2, dad1, dad2, kid1,kid2, path, summary_name;
	ifstream vcf, dad_kid_file;
	ofstream consensus, output, summary, mom;
	vector<comparisons> dad_kid;
	vector<locus_data> locus_info;
	bool first = true;

	path = "../../../results/biallelic/";
	dad_kid_name = "../../../results/dad.kid.pairs.fullnames.txt";
	vcf_name = "../../../results/stacks/batch_1.vcf";
	summary_name = path + "biallelic_maternal.vcf";

	
	dad_kid_file.open(dad_kid_name);
	FileTest(dad_kid_file, dad_kid_name);
	index = 0;
	while (universal_getline(dad_kid_file, line))
	{
		if (!dad_kid_file.eof())
		{
			dad_kid.push_back(comparisons());
			stringstream ss;
			ss.str(line);
			ss >> dad_kid[index].dad_name;
			ss >> dad_kid[index].kid_name;
			index++;
		}
	}
	dad_kid_file.close();

	cout << "\nThere are " << dad_kid.size() << " dad-kid comparisons.\n";
	
	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	summary.open(summary_name);
	summary << "##source='Includes Inferred Maternal Alleles'";
	locus_index = line_count = snp_count = 0;
	while (universal_getline(vcf, line))
	{
		if (!vcf.eof())
		{
			if (line.substr(0, 2) == "##")
				summary << '\n' << line;
			if (line.substr(0, 2) == "#C")
			{
				//need to determine male-offspring indices
				stringstream ss;
				ss.str(line);
				count = index = 0;
				summary << '\n';
				while (ss >> stemp)
				{
					if (count > 8)
					{
						for (i = 0; i < dad_kid.size(); i++)
						{
							if (dad_kid[i].dad_name == stemp)
								dad_kid[i].dad_index = index;
							if (dad_kid[i].kid_name == stemp)
								dad_kid[i].kid_index = index;
						}
						index++;
					}
					else
					{
						if (count == 0)
							summary << stemp;
						else
							summary << '\t' << stemp;
					}
					count++;
				}//while ss >>stemp
				for (i = 0; i < dad_kid.size(); i++)
				{
					summary << '\t' << dad_kid[i].dad_name << '\t' << dad_kid[i].kid_name << '\t' << "MOM" << dad_kid[i].dad_name.substr(3);
				}
			}
			if (line.substr(0, 1) != "#") //then it's a new locus
			{

				stringstream ss;
				ss.str(line);
				ss >> chrom >> pos >> ID >> ref >> alt >> quality >> filter >> info >> format;
				locus_info.push_back(locus_data());
				locus_info[locus_index].pos = pos;
				locus_info[locus_index].chrom = chrom;
				locus_info[locus_index].ref = ref;
				locus_info[locus_index].alt = alt;
				locus_info[locus_index].id = atoi(ID.c_str());
				snp_count++;
				summary << '\n' << chrom << '\t' << pos << '\t' << ID << '\t' << ref << '\t' << alt << '\t' << quality << '\t' << filter << '\t' << info << '\t' << "GT";
				stringstream ssi;
				ssi << format;
				count = 0;
				vector<string> format_fields;
				while (getline(ssi, temp, ':'))
				{
					if (temp == "GT") //determine which slot in the FORMAT is the genotype (GT)
						gt_index = count;
					count++;
				}
				
				index = 0;
				vector<string> alleles;
				while (ss >> stemp)
				{
					stringstream ssii;
					ssii << stemp;
					vector<string> vectemp;
					while (getline(ssii, temp, ':'))
					{
						vectemp.push_back(temp);
					}
					alleles.push_back(vectemp[gt_index]);
				}
				//infer maternal alleles
				for (i = 0; i < dad_kid.size(); i++)
				{
					stringstream ssiii;
					ssiii.str(alleles[dad_kid[i].dad_index]);
					vector<string> vectemp;
					while (getline(ssiii, temp, '/'))
					{
						vectemp.push_back(temp);
					}
					dad1 = vectemp[0];
					dad2 = vectemp[1];
					stringstream ssiv;
					ssiv.str(alleles[dad_kid[i].kid_index]);
					vectemp.resize(0);
					while (getline(ssiv, temp, '/'))
					{
						vectemp.push_back(temp);
					}
					kid1 = vectemp[0];
					kid2 = vectemp[1];
					
					string mom_allele = ".";
					if (dad1 == dad2 && kid1 == kid2) //the case where both are homozygous
					{
						if (dad1 == kid1)
							mom_allele = dad1;
					}
					if (dad1 == dad2 && kid1 != kid2)//the case where dad is homozygous but kid is het
					{
						if (dad1 == kid1)
							mom_allele = kid2;
						if (dad1 == kid2)
							mom_allele = kid1;
					}
					if (dad1 != dad2 && kid1 == kid2) //the case where dad is het but off is hom
					{
						if (dad1 == kid1 || dad2 == kid1)
							mom_allele = kid1;
					}
					if (dad1 != dad2 && kid1 != kid2)//if they're both hets you can't do anything with it.
						mom_allele = ".";
					if (dad1 == "." || kid1 == ".")
						mom_allele = ".";
					dad_kid[i].mom_allele.push_back(mom_allele);
					summary << '\t' << dad1 << "/" << dad2 << '\t' << kid1 << "/" << kid2 << '\t' << mom_allele << "/" << mom_allele;
			
				}	//dad_kid.size()
				if (snp_count == 20000)
				{
					cout << "\n20000 SNP entries finished. Writing those to files.\n";
					if (first == true)
					{
						for (i = 0; i < dad_kid.size(); i++)
						{
							string mom_name;
							if (dad_kid[i].kid_name.size() > 16)
							{
								int second_und = dad_kid[i].kid_name.size() - 16;
								mom_name = path + "MOM" + dad_kid[i].kid_name.substr(10, second_und) + ".bi.txt";
							}
							else
							{
								mom_name = path + "MOM" + dad_kid[i].kid_name.substr(3, dad_kid[i].kid_name.size()) + ".bi.txt";
							}
							mom.open(mom_name);
							mom << "SNPID\tAllele1\tAllele2";
							for (ii = 0; ii < locus_info.size(); ii++)
							{
								if (dad_kid[i].mom_allele[ii] != ".")//don't write empty ones to file
								{
									mom << '\n' << locus_info[ii].id << "." << locus_info[ii].pos << '\t';
									if (dad_kid[i].mom_allele[ii] == "0")
										mom << locus_info[ii].ref << "\t-";
									else
										mom << locus_info[ii].alt << "\t-";
								}
							}
							mom.close();
						}
					}
					else
					{
						for (i = 0; i < dad_kid.size(); i++)
						{
							string mom_name;
							if (dad_kid[i].kid_name.size() > 16)
							{
								int second_und = dad_kid[i].kid_name.size() - 16;
								mom_name = path + "MOM" + dad_kid[i].kid_name.substr(10, second_und) + ".bi.txt";
							}
							else
							{
								mom_name = path + "MOM" + dad_kid[i].kid_name.substr(3, dad_kid[i].kid_name.size()) + ".bi.txt";
							}
							mom.open(mom_name, ios::app);
							for (ii = 0; ii < locus_info.size(); ii++)
							{
								if (dad_kid[i].mom_allele[ii] != ".")//don't write empty ones to file
								{
									mom << '\n' << locus_info[ii].id << "." << locus_info[ii].pos << '\t';
									if (dad_kid[i].mom_allele[ii] == "0")
										mom << locus_info[ii].ref << "\t-";
									else
										mom << locus_info[ii].alt << "\t-";
								}
							}
							mom.close();
							dad_kid[i].mom_allele.resize(0);
						}
					}
					first = false;
					snp_count = 0;
					locus_index = -1;
					locus_info.resize(0);
				}//snp_count = 20000
				line_count++;
				locus_index++;
				if (line_count % 1000 == 0)
					cout << "\nCompleted " << line_count << " SNP entries.";
			}//if #		
		}//eof
	}//getline
	vcf.close();
	summary.close();

	for (i = 0; i < dad_kid.size(); i++)
	{
		string mom_name;
		if (dad_kid[i].kid_name.size() > 16)
		{
			int second_und = dad_kid[i].kid_name.size() - 16;
			mom_name = path + "MOM" + dad_kid[i].kid_name.substr(10, second_und) + ".bi.txt";
		}
		else
		{
			mom_name = path + "MOM" + dad_kid[i].kid_name.substr(3, dad_kid[i].kid_name.size()) + ".bi.txt";
		}
		mom.open(mom_name, ios::app);
		for (ii = 0; ii < locus_info.size(); ii++)
		{
			if (dad_kid[i].mom_allele[ii] != ".")//don't write empty ones to file
			{
				mom << '\n' << locus_info[ii].id << "." << locus_info[ii].pos << '\t';
				if(dad_kid[i].mom_allele[ii] == "0")
					mom << locus_info[ii].ref << "\t-";
				else
					mom << locus_info[ii].alt << "\t-";
			}
		}
		mom.close();
		dad_kid[i].mom_allele.resize(0);
	}
	
	cout << "\nSuccessfully completed " << line_count << " SNPs.";
		
	cout << "\nInput integer to quit.\n";
	cin >> i;
	return 0;
}