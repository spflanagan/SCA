//Author: Sarah P. Flanagan
//Date: 10 December 2015
//filter_vcf
//this program takes GATK vcf after variant call pipeline, minus filtering step
//it will filter following Monnahan et al. (2015)
//applies following filters:
//1. Keeps SNPs with two segregating alleles
//2. Present in at least 50 parents and 50 progeny
//3. initial estimated allele frequency between 0.05 and 0.95
//4. GATK haplotype score < 13
//5. Mapping quality score >= 30
//6. Average read coverage per individual >= 1 and <= 100

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
	int i, pos, line_count, count, prog_start, prog_end, prog_count, adult_count, ind_count;
	string vcf_name, vcf_out_name, line, summary_name;
	string chrom, ID, ref, alt, filter, info, format, temp, temp2, last_chr, quality;
	double avg_dp, avg_gq;//quality, 
	ifstream vcf;
	ofstream vcf_out, summary;
	bool pass_initial_screen, pass_second_screen;
	int dpi, gqi, gti, keep_count;

	vcf_name = "E://ubuntushare//SCA//batch_1.vcf";
	vcf_out_name = "batch_1_filtered.vcf";
	summary_name = "batch1_filtered_summary.txt";

	vcf.open(vcf_name);
	FileTest(vcf, vcf_name);
	vcf_out.open(vcf_out_name);
	summary.open(summary_name);
	summary << "Scaffold\tPos\tAvgDP\tAvgGQ\tProgCount\AdultCount";
	line_count = keep_count = 0;
	while (universal_getline(vcf, line))
	{
		if (!vcf.eof())
		{
			if (line.substr(0, 2) != "##")
			{
				stringstream ss;
				ss.str(line);
				if (line.substr(0, 2) == "#C")
				{//evaluate individual indices for adults and offspring
					vcf_out << '\n' << line;
					count = 0;
					bool prog_found = false;
					prog_start = prog_end = -5;
					while (ss >> temp)
					{
						if (temp.substr(0, 3) == "OFF")
							prog_found = true;
						else 
							prog_found = false;
						if (prog_found == true && prog_start < 0)
							prog_start = count;
						if (prog_found == false && prog_start > 0 && prog_end < 0)
							prog_end = count - 1;
						count++;
					}
				}
				else
				{//then it's a real line of info
					ss >> chrom >> pos >> ID >> ref >> alt >> quality >> filter >> info >> format;
					//break down info
					vector<string> inds;
					stringstream ssa;
					ssa << alt;
					vector<string> alt_temp;
					while (getline(ssa, temp, ','))
						alt_temp.push_back(temp);
					if (alt_temp.size() < 2)
					{
						stringstream ssi;
						ssi << info;
						pass_initial_screen = true;
						while (getline(ssi, temp, ';'))
						{//check depth in individuals
							if (temp.substr(0, 2) == "AF")
							{//Allele frequency between 0.05 and 0.95
								if (atof(temp.substr(3, temp.size()).c_str()) <= 0.05 || atof(temp.substr(3, temp.size()).c_str()) >= 0.95)
									pass_initial_screen = false;
							}
							if (temp.substr(0, 2) == "MQ")
							{//must have mapping quality >=30
								if (atof(temp.substr(3, temp.size()).c_str()) <= 30)
									pass_initial_screen = false;
							}
						}
						if (pass_initial_screen == true)
						{
							stringstream ssi;
							ssi << format;
							dpi = gqi = gti = count = 0;
							while (getline(ssi, temp, ':'))
							{
								if (temp == "DP")
									dpi = count;
								if (temp == "GQ")
									gqi = count;
								if (temp == "GT")
									gti = count;
								count++;
							}
							avg_dp = avg_gq = 0;
							ind_count = prog_count = adult_count = 0;
							while (ss >> temp)
							{
								inds.push_back(temp);
								stringstream ssii;
								ssii << temp;
								count = 0;
								while (getline(ssii, temp2, ':'))
								{
									if (count == dpi)
										avg_dp = avg_dp + atof(temp2.c_str());
									if (count == gqi)
										avg_gq = avg_gq + atof(temp2.c_str());
									if (count == gti)
									{
										if (temp2.substr(0, 1) != "." || temp2.substr(1, 1) != ".")
										{
											if (ind_count >= prog_start && ind_count <= prog_end)
												prog_count++;
											else
												adult_count++;
										}
									}
									count++;
								}//while reading individual
								ind_count++;
							}//while reading all individuals
							avg_dp = avg_dp / ind_count;
							avg_gq = avg_gq / ind_count;
							pass_second_screen = true;
							if (avg_dp < 1 || avg_dp > 100)
								pass_second_screen = false;
							if (avg_gq >= 13)
								pass_second_screen = false;
							if (prog_count < 50 && adult_count < 50)
								pass_second_screen = false;
							if (pass_second_screen == true)
							{
								//if it passes all filters, send it to the new file!
								vcf_out << '\n' << line;
								keep_count++;
								if (keep_count % 1000 == 0)
									cout << keep_count << " SNPs retained\n";
								summary << '\n' << chrom << '\t' << pos << '\t' << avg_dp << '\t' << avg_gq << '\t' << prog_count << '\t' << adult_count;
							}
						}
					}
				}
			}//if ##
			else
			{
				if (line_count == 0)
					vcf_out << line;
				else
					vcf_out << '\n' << line;
			}
			line_count++;
		}
	}//getline
	vcf.close();
	vcf_out.close();
	summary.close();

	cout << "A total of " << keep_count << " SNPs passed the filtering steps.\n";

	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}

