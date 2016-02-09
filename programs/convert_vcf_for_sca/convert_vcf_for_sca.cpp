//Author: Sarah P. Flanagan
//Date: 15 December 2015
//Purpose: This code will take a vcf file (preferably pre-filtered) and convert it 
//The output format is compatible with Monnahan et al. (2015)'s selection components analysis python script
//(I think)
//It has one row per family group, all individuals grouped per SNP
//there can be different numbers of columns per row.


#include <iostream>
#include <fstream>
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


int main()
{
	int i, ii, line_count, snp_count, count, prog_start, prog_end;
	string vcf_name, output_name, line, stemp;
	string chrom, ID, ref, alt, filter, info, format, temp, temp2, last_chr;
	int pos, index, num_removed = 0;
	string quality;
	ifstream vcf;
	ofstream output;
	vector <string> keep_format;
	keep_format = { "GT", "AD", "DP", "GQ", "PL" };

	vcf_name = "E://ubuntushare//SCA//monnahan//genotype_output_filtered.recode.vcf";
	output_name = "E://ubuntushare//SCA//monnahan//genotype_output_filtered.vcf";

	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	output.open(output_name);
	line_count = snp_count = 0;
	while (universal_getline(vcf, line))
	{
		if (!vcf.eof())
		{
			if (line.substr(0, 1) != "#")
			{
				stringstream ss;
				ss.str(line);
				ss >> chrom >> pos >> ID >> ref >> alt >> quality >> filter >> info >> format;
				output << '\n' << chrom << '\t' << pos << '\t' << ID << '\t' << ref << '\t' << alt << '\t' << quality << '\t' << filter << '\t' << info << '\t';
				output << "GT:AD:DP:GQ:PL";
				stringstream ssi;
				ssi << format;
				index = count = 0;
				vector<string> format_fields;
				bool fix_format = false;
				while (getline(ssi, temp, ':'))
					format_fields.push_back(temp);
				if (format_fields.size() < 5)
					fix_format = true;
					
				if (fix_format == true)
				{
					count = 0;
					while (ss >> stemp)
					{
						stringstream ssii;
						ssii << stemp;
						count = index = 0;
						output << '\t';
						vector<string> vectemp;
						while (getline(ssii, temp, ':'))
						{
							vectemp.push_back(temp);
						}
						if (vectemp[0] == "./.")
							output << "./.";
						else
						{
							for (i = 0; i < keep_format.size(); i++)
							{
								index = -1;
								for (ii = 0; ii < format_fields.size(); ii++)
								{
									if (keep_format[i] == format_fields[ii])
									{
										index = ii;
									}
								}
								if (index >= 0)
								{
									if (keep_format[i] == "PL")
									{
										vector<string> gt;
										stringstream ssiii;
										ssiii << vectemp[index];
										while (getline(ssiii, temp, ','))
											gt.push_back(temp);
										for (int j = 0; j < gt.size(); j++)
										{
											if (gt[j] == ".")
												output << "0";
											else
												output << gt[j];
											if (j < gt.size() - 1)
												output << ",";
										}
										
									}
									else
									{
										output << vectemp[index] << ":";
										index++;
									}
								}
								else
								{
									if (keep_format[i] == "GT")
										output << "./.:";
									if (keep_format[i] == "AD")
										output << ".,.:";
									if (keep_format[i] == "DP")
										output << ".:";
									if (keep_format[i] == "GQ")
										output << ".:";
									if (keep_format[i] == "PL")
										output << ".,.,.";									
								}
							}
						}//end else
						count++;
					}//end while ss >> temp
					num_removed++;
				}
				else
				{
					while (ss >> stemp)
					{
						stringstream ssii;
						ssii << stemp;
						count = index = 0;
						output << '\t';
						vector<string> vectemp;
						while (getline(ssii, temp, ':'))
						{
							vectemp.push_back(temp);
						}
						if (vectemp[0] == "./.")
							output << "./.";
						else
						{
							for (i = 0; i < vectemp.size() - 1; i++)
							{
								output << vectemp[i] << ":";
							}
							output << vectemp[vectemp.size() - 1];
						}
					}
				}
			}//if #
			else
			{
				if (line_count > 0)
					output << '\n';
				output << line;
			}
			line_count++;

			if (line_count % 1000 == 0)
				cout << "\nCompleted " << line_count << " SNP entries.";
		}
	}//getline
	vcf.close();
	output.close();
	cout << num_removed << " loci were corrected.\n";
	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}