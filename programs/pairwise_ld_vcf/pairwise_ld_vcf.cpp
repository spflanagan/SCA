//Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
//Date: 9 May 2016
//Last updated: 9 May 2016
//Purpose: To calculate pairwise LD between all loci on a chromosome, using a vcf file as input.

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
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

class ld_info
{
public: 
	double d, dprime;
	string chrom1, chrom2, locus1, locus2;
	double pA, pB, pAB, paB, pAb, pab;

	ld_info()
	{
		d = dprime = double();
		chrom1 = chrom2 = locus1 = locus2 = string();
		pA = pB = pAB = paB = pAb = pab = double();
	}
};

class locus_info
{
public:
	string chrom, locusID;
	int bp;
	vector<int> maternal;
	vector<int> paternal;

	locus_info()
	{
		chrom = locusID = string();
		bp = int();
		maternal = paternal = vector<int>();
	}
};

ld_info AdultPopLD(locus_info & locus1, locus_info &locus2)
{
	//To do this, need to compare allele frequencies. 
	//D=x11-p1q1
	//x11 = observed freq of major allele in locus A and major allele in locus B
	//p1 = observed freq of major allele in locus A
	//q1 = observed freq of major allele in locus B
	//if D > 0, Dmax = min(p1q1, p2q2)
	//if D < 0, Dmax = min(p1q2, p2q1)
	//D' = D/Dmax
	ld_info result;
	int NumAlleles = 2;
	int num_adults, nad1, nad2;
	vector <int> num_allele_1;
	vector <int> num_allele_2;
	vector <double> freq_allele_1;
	vector <double> freq_allele_2;
	double joint_alleles[2][2];
	double D[2][2]; 
	double Dmax[2][2];
	//double **joint_alleles = new double *[NumAlleles];
	//double **D = new double *[NumAlleles];
	//double **Dmax = new double *[NumAlleles];
	int count_a = 0;
	int f, ff, fff, al, count;
	double Dprime, d_allele_avgs;
	int numA1B1 = 0;

	num_adults = nad1 = nad2 = 0;
	/*for (al = 0; al < NumAlleles; al++)
	{
		joint_alleles[al] = new double[NumAlleles];
		D[al] = new double[NumAlleles];
		Dmax[al] = new double[NumAlleles];
	}*/
	//figure out which loci are polymorphic
	for (al = 0; al < NumAlleles; al++)
	{
		num_allele_1.push_back(0);
		num_allele_2.push_back(0);
		freq_allele_1.push_back(0);
		freq_allele_2.push_back(0);
		for (f = 0; f < NumAlleles; f++)
		{
			joint_alleles[al][f] = 0;
			D[al][f] = 0;
			Dmax[al][f] = 0;
		}
	}

	int maternal_1, maternal_2, paternal_1, paternal_2;

	for (fff = 0; fff < locus1.maternal.size(); fff++)
	{
		maternal_1 = locus1.maternal[fff];
		maternal_2 = locus2.maternal[fff];
		paternal_1 = locus1.paternal[fff];
		paternal_2 = locus2.paternal[fff];
		/*if (maternal_1 >= 0)
		{
			num_allele_1[maternal_1]++;
			num_allele_1[paternal_1]++;
			nad1++;
		}
		if (paternal_2 >= 0)
		{
			num_allele_2[maternal_2]++;
			num_allele_2[paternal_2]++;
			nad2++;
		}*/
		if (maternal_1 >= 0 && maternal_2 >= 0)
		{			
			num_allele_1[maternal_1]++;
			num_allele_1[paternal_1]++;
			nad1++;
			num_allele_2[maternal_2]++;
			num_allele_2[paternal_2]++;
			nad2++; 
			joint_alleles[maternal_1][maternal_2]++;
			joint_alleles[paternal_1][paternal_2]++;
			num_adults++;
		}
		
	}//end for fff < PopulationSize

	int major_allele_1 = 0;
	int major_allele_2 = 0;
	for (al = 0; al < NumAlleles; al++)
	{
		freq_allele_1[al] = (num_allele_1[al]) / (2 * double(nad1));
		freq_allele_2[al] = (num_allele_2[al]) / (2 * double(nad2));
		if (freq_allele_1[al] > freq_allele_1[major_allele_1])
			major_allele_1 = al;
		if (freq_allele_2[al] > freq_allele_2[major_allele_2])
			major_allele_2 = al;
		for (f = 0; f < NumAlleles; f++)
			joint_alleles[al][f] = joint_alleles[al][f] / (2 * double(num_adults));
	}

	if (freq_allele_1[major_allele_1] != 1 && freq_allele_2[major_allele_1] != 1)//if it's polymorphic
	{
		d_allele_avgs = 0;
		count = 0;
		for (f = 0; f < NumAlleles; f++)
		{
			for (ff = 0; ff < NumAlleles; ff++)
			{
				if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
				{
					D[f][ff] = joint_alleles[f][ff] - freq_allele_1[f] * freq_allele_2[ff];
					d_allele_avgs = d_allele_avgs + fabs(D[f][ff]);
					count++;
					if (D[f][ff] < 0)
						Dmax[f][ff] = min(freq_allele_1[f] * freq_allele_2[ff], (1 - freq_allele_1[f])*(1 - freq_allele_2[ff]));
					else
						Dmax[f][ff] = min((1 - freq_allele_1[f])*freq_allele_2[ff], freq_allele_1[f] * (1 - freq_allele_2[ff]));
				}
			}
		}
		double dcount = count;
		d_allele_avgs = d_allele_avgs / dcount;

		Dprime = 0;
		bool decentDmax = true;
		for (f = 0; f < NumAlleles; f++)
		{
			for (ff = 0; ff < NumAlleles; ff++)
			{
				if (freq_allele_1[f] > 0 && freq_allele_2[ff] > 0)
				{
					if (Dmax[f][ff] > 0)
						Dprime = Dprime + freq_allele_1[f] * freq_allele_2[ff] * fabs(D[f][ff]) / Dmax[f][ff];
					else
						decentDmax = false;

				}
			}
		}
		if (!decentDmax)
			Dprime = -5;

	}

	if (Dprime > 1.1)
		cout << "D' greater than 1.";
	result.d = d_allele_avgs;
	result.dprime = Dprime;
	result.chrom1 = locus1.chrom;
	result.locus1 = locus1.locusID;
	result.chrom2 = locus2.chrom;
	result.locus2 = locus2.locusID;
	
	//pA=freq(0 loc 1), pB = freq(0 loc 2)
	result.pA = freq_allele_1[0];
	result.pB = freq_allele_2[0];
	result.pAB = joint_alleles[0][0];
	result.pAb = joint_alleles[0][1];
	result.paB = joint_alleles[1][0];
	result.pab = joint_alleles[0][0];

	/*for (f = 0; f < NumAlleles; f++)
	{
		delete[] joint_alleles[f];
		delete[] D[f];
		delete[] Dmax[f];
	}
	delete[] joint_alleles;
	delete[] D;
	delete[] Dmax;*/

	return result;
}//end Adult Pop LD

int main()
{
	int i, ii, iii, count;
	string vcf_name, out_name,matrix_name, any_matrix_name,line, temp, last_lg, this_lg;
	ifstream vcf;
	ofstream out, matrix, any_matrix;
	vector<locus_info> vcf_dat;
//	vector<ld_info> ld;

	vcf_name = "../../results/stacks/batch_1.adults_MAF1.recode.vcf";
	out_name = "../../results/genome_paper/ld_info_adults_maf1.txt";
	matrix_name = "../../results/genome_paper/ld_matrix_adults_maf1_lg1.txt";
	string any_matrix_base = "../../results/genome_paper/ld_matrix_adults_maf1_";

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
				ss >> vcf_dat[count].chrom >> vcf_dat[count].locusID >> vcf_dat[count].bp >> temp >> temp >> temp >> temp >> temp >> temp;
				while (ss >> temp)
				{
					if (temp.substr(0,3) != "./.")
					{
						string allele1 = temp.substr(0, 1);
						string allele2 = temp.substr(2, 1);
						vcf_dat[count].maternal.push_back(atoi(allele1.c_str()));
						vcf_dat[count].paternal.push_back(atoi(allele2.c_str()));
					}
					else
					{
						vcf_dat[count].maternal.push_back(-5);
						vcf_dat[count].paternal.push_back(-5);
					}
				}
				count++;
			}//if it's a data-filled line
		}//while vcf line
	}//while vcf
	vcf.close();
	cout << "Found " << count << " locus records.\n";

	out.open(out_name);
	out << "Chrom\tLocIDA\tpA\tLocIDB\tpB\tpAB\tpAb\tpaB\tpab\tD\tD'";
	matrix.open(matrix_name);
	this_lg = vcf_dat[0].chrom;
	cout << this_lg;
	any_matrix_name = any_matrix_base + this_lg + ".txt";
	any_matrix.open(any_matrix_name);
	for (i = 0; i < vcf_dat.size(); i++)
	{
		if (vcf_dat[i].chrom == "LG1")
			matrix << '\t' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		if (vcf_dat[i].chrom == this_lg)
			any_matrix << '\t' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		else
		{
			any_matrix.close();
				last_lg = this_lg;
				this_lg = vcf_dat[i].chrom;
				//cout << this_lg << '\n';
				any_matrix_name = any_matrix_base + this_lg + ".txt";
				any_matrix.open(any_matrix_name, std::ios_base::app);
				any_matrix << '\t' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		}
	}
	count = 0;
	this_lg = vcf_dat[0].chrom;
	any_matrix_name = any_matrix_base + this_lg + ".txt";
	any_matrix.open(any_matrix_name);
	for (i = 0; i < vcf_dat.size(); i++)
	{
		if (vcf_dat[i].chrom == "LG1")
			matrix << '\n' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		if (vcf_dat[i].chrom == this_lg)
			any_matrix << '\n' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		else
		{
			any_matrix.close();
				last_lg = this_lg;
				this_lg = vcf_dat[i].chrom;
				any_matrix_name = any_matrix_base + this_lg + ".txt";
				any_matrix.open(any_matrix_name, std::ios_base::app);
				any_matrix << '\n' << vcf_dat[i].chrom << "." << vcf_dat[i].locusID << "." << vcf_dat[i].bp;
		}
		for (ii = 0; ii < vcf_dat.size(); ii++)
		{
			
			if (i != ii && vcf_dat[i].chrom == vcf_dat[ii].chrom)
			{
				ld_info ld;
				ld = AdultPopLD(vcf_dat[i], vcf_dat[ii]);
				out << '\n' << ld.chrom1 << '\t' << ld.locus1 << '\t' << ld.pA << '\t'
					<< ld.locus2 << '\t' << ld.pB << '\t' << ld.pAB << '\t'
					<< ld.pAb << '\t' << ld.paB << '\t' << ld.pab << '\t'
					<< ld.d << '\t' << ld.dprime;
				
				if (vcf_dat[i].chrom == "LG1" && vcf_dat[ii].chrom == "LG1")
					matrix << '\t' << ld.dprime;
				if (vcf_dat[i].chrom == this_lg && vcf_dat[ii].chrom == this_lg)
					any_matrix << '\t' << ld.dprime;
				
				count++;
			}
			else
			{
				if (vcf_dat[i].chrom == "LG1" && vcf_dat[ii].chrom == "LG1")
					matrix << '\t' << -1;
				if (vcf_dat[i].chrom == this_lg && vcf_dat[ii].chrom == this_lg)
					any_matrix << '\t' << -1;
				
			}
		}
	}
	out.close();
	matrix.close();
	cout << "\nPerformed " << count << " pairwise LD calculations.";
	cout << "\nDone! Input integer to quit\n";
	cin >> i;
	return 0;
}