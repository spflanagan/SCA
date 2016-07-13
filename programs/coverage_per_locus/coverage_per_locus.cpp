//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date: 16 June 2016
//Purpose: Calculate the mean and variance in coverage per locus using matches files from Stacks


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

class locus
{
public:
	int cat_id;
	int count;
	vector<string> alleles;
	vector<int> allele_count;

	locus()
	{
		count = int();
		cat_id = int();
		alleles = vector < string >();
		allele_count = vector<int>();
	}
};

class locus_info
{
public:
	string chrom, locusID;
	int bp;
	vector<char> maternal;
	vector<char> paternal;
	vector<int> depth1;
	vector<int> depth2;

	locus_info()
	{
		chrom = locusID = string();
		bp = int();
		maternal = paternal = vector<char>();
		depth1 = depth2 = vector<int>();
	}
};

int main()
{
	int i, ii, iii, count,het_count;
	char alleles[2];
	string vcf_name, line, temp,t2,tt2, out_name;
	ifstream vcf, popmap;
	ofstream output;
	vector<locus_info> vcf_dat;
	vector<string> vcf_ind_ids;
	vector<string> pops;
	vector<int> locus_depth;
	vector<double> prop_het;
	vector<double> het_depth_ratio;

	vcf_name = "../../results/stacks/batch_1.vcf";
	out_name = "../../results/quality_control/coverage_info.txt";

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
				while (ss >> temp)
				{
					if (temp.substr(0, 3) != "./.")
					{
						string allele1 = temp.substr(0, 1);
						string allele2 = temp.substr(2, 1);
						stringstream t;
						t.str(temp);
						for(i = 0; i < 3;i++)//the third one is the AD
							getline(t, t2,':');
						stringstream tt;
						tt.str(t2);
						getline(tt, tt2, ',');
						int d1 = stoi(tt2);
						getline(tt, tt2, ',');
						int d2 = stoi(tt2);
						char al1 = alleles[atoi(allele1.c_str())];
						char al2 = alleles[atoi(allele2.c_str())];
						vcf_dat[count].maternal.push_back(al1);
						vcf_dat[count].paternal.push_back(al2);
						vcf_dat[count].depth1.push_back(d1);
						vcf_dat[count].depth2.push_back(d2);
					}
					else
					{
						vcf_dat[count].maternal.push_back(0);
						vcf_dat[count].paternal.push_back(0);
						vcf_dat[count].depth1.push_back(0);
						vcf_dat[count].depth2.push_back(0);
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
	cout << "Found " << count << " locus records.\n";

	//now get coverage per locus
	output.open(out_name);
	cout << "\nWriting to " << out_name;
	output << "Chrom\tBP\tLocus\tTotalDepth\tAvgDepthPerInd\tVarianceInDepthPerInd\tPropHet\tHetDepthRatio\tNumHet\tTotalN";
	for (i = 0; i < vcf_dat.size(); i++)
	{
		locus_depth.push_back(0);
		prop_het.push_back(0);
		het_depth_ratio.push_back(0);
		double depth_variance = 0;
		double depth_mean = 0;
		count = het_count= 0;
		for (ii = 0; ii < vcf_dat[i].depth1.size(); ii++)
		{
			locus_depth[i] = locus_depth[i] + vcf_dat[i].depth1[ii] + vcf_dat[i].depth2[ii];
			if (vcf_dat[i].maternal[ii] != vcf_dat[i].paternal[ii]) //they're heterozygotes
			{
				//prop_het[i]++;
				double d1 = vcf_dat[i].depth1[ii];
				double d2 = vcf_dat[i].depth2[ii];
				het_depth_ratio[i] = het_depth_ratio[i] + (d1 /d2);
				het_count++;
			}
			if (vcf_dat[i].maternal[ii] != '\0')
				count++;
		}
		prop_het[i] = double(het_count) / double(count);
		if (het_count >0)
			het_depth_ratio[i] = double(het_depth_ratio[i]) / double(het_count);
		else
			het_depth_ratio[i] = 0;
		//calculate variance in depth
		depth_mean = locus_depth[i] / double(count);
		for (ii = 0; ii < vcf_dat[i].depth1.size(); ii++)
		{
			if (vcf_dat[i].maternal[ii] != '\0')
			{
				int depth = vcf_dat[i].depth1[ii] + vcf_dat[i].depth2[ii];
				depth_variance = depth_variance + (depth - depth_mean)*(depth - depth_mean);
				count++;
			}
		}
		depth_variance = depth_variance / double(count);
		output << '\n' << vcf_dat[i].chrom << '\t' << vcf_dat[i].bp << '\t' << vcf_dat[i].locusID << '\t' <<
			locus_depth[i] << '\t' << depth_mean << '\t' << depth_variance << '\t' << prop_het[i] << '\t' << het_depth_ratio[i] <<
			'\t' << het_count << '\t' << count;
	}
	output.close();

	

	cout << "\nDone!";
	cin >> i;
	return 0;
}