//Author: Sarah P. Flanagan
//Date: 21 September 2015
//Purpose: Use allele info (from process_alleles_files or otherwise) and add infer maternal genotypes from paired dads and offspring

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

class summary_structure
{
public:
	vector <int> locus_id;
	vector <vector<string>> all_alleles;
	vector <int> num_in_dads;
	vector <int> num_matched;

	summary_structure()
	{
		locus_id = vector<int>();
		all_alleles = vector<vector<string>>();
		num_in_dads = vector<int>();
		num_matched = vector<int>();
	}
};

class catalog
{
public:
	int catalog_id;
	vector<string> alleles;

	catalog()
	{
		catalog_id = int();
		alleles = vector<string>();
	}
};

class alleles
{
public:
	string hap;
	int count;

	alleles()
	{
		hap = string();
		count = int();
	}
};

class genotype
{
public:
	vector <int> locus;
	vector <int> catalog;
	vector <alleles> hap1;
	vector <alleles> hap2;

	genotype()
	{
		catalog = vector<int>();
		locus = vector<int>();
		hap1 = vector<alleles>();
		hap2 = vector<alleles>();
	}
};

class one_pair
{
public:
	string dad_id, off_id;
	genotype dad_alleles;
	genotype off_alleles;
	vector<string> mom_allele;

	one_pair()
	{
		dad_id = off_id = string();
		dad_alleles = genotype();
		off_alleles = genotype();
		mom_allele = vector<string>();
	}
};

class one_match
{
public:
	vector<string> haplotype;
	int cat_id, loc_id;
	vector<int> count;

	one_match()
	{
		haplotype = vector<string>();
		cat_id = int();
		loc_id = int();
		count = vector<int>();
	}
};

int main()
{
	int end, count, allele_count, line_count;
	int tmp1, tmp2, id, samp_id, loc_id, depth, last_allele;
	string haplotype;
	double log_like, perc;
	string paired_list_name, dad_name, off_name, in_path, out_path, consensus_name, output_name, summary_name, mom_name, cat_alleles_name;
	string line, line1, line2;
	ifstream paired_list, dad, off, cat_alleles;
	ofstream consensus, output, summary, mom;
	string a, b;
	int i, ii;
	summary_structure consensus_loci;
	vector<catalog> cat;

	in_path = "../../results/stacks/"; 
	out_path = "../../results/haplotypes/";
	cat_alleles_name = in_path + "batch_1.catalog.alleles.tsv";
	paired_list_name = out_path + "dad.kid.pairs.txt";
	summary_name = out_path + "infer_maternal_summary.txt";

	////read in catalog
	//cat_alleles.open(cat_alleles_name);
	//FileTest(cat_alleles, cat_alleles_name);
	//last_allele = 0;
	//while (universal_getline(cat_alleles, line))
	//{
	//	if (!cat_alleles.eof())
	//	{
	//		stringstream ss;
	//		ss.str(line);
	//		ss >> tmp1 >> tmp2 >> id >> haplotype >> perc >> count;
	//		if (id != last_allele)//assumes the ids are in order
	//		{
	//			cat.push_back(catalog());
	//			cat.back().catalog_id = id;
	//			cat.back().alleles.push_back(haplotype);
	//		}
	//		else
	//			cat.back().alleles.push_back(haplotype);
	//	}
	//}
	//cat_alleles.close();

	summary.open(summary_name);
	summary << "Parent1ID\tOffspringID\tNumLociParent1\tNumInferredLoci";
	paired_list.open(paired_list_name);
	FileTest(paired_list, paired_list_name);
	ii = 0;
	while (universal_getline(paired_list, line))//each pair
	{
		if (!paired_list.eof())
		{
			stringstream ss;
			ss.str(line);
			ss >> dad_name >> off_name;
			one_pair dad_off;
			dad_off.dad_id = dad_name;
			dad_off.off_id = off_name;
			summary << '\n' << dad_name << '\t' << off_name;
			dad_name = in_path + "sample_" + dad_name + "_align.matches.tsv";
			off_name = in_path + "sample_" + off_name + "_align.matches.tsv";
			dad.open(dad_name);
			FileTest(dad, dad_name);
			line_count = count = 0;
			vector<one_match> temp;
			while (universal_getline(dad, line1))
			{
				if (line_count == 0)
				{
					stringstream ss1;
					ss1.str(line1);
					ss1 >> tmp1 >> tmp2 >> id >> samp_id >> loc_id >> haplotype >> depth >> log_like;
					temp.push_back(one_match());
					temp[0].cat_id = id;
					temp[0].count.push_back(depth);
					temp[0].haplotype.push_back(haplotype);
					temp[0].loc_id = loc_id;
				}
				else
				{
					stringstream ss1;
					ss1.str(line1);
					ss1 >> tmp1 >> tmp2 >> id >> samp_id >> loc_id >> haplotype >> depth >> log_like;
					allele_count = 0;
					bool found = false;
					for (size_t t = 0; t < temp.size(); t++)
					{
						if (temp[t].cat_id == id)
						{
							found = true;
							temp[t].haplotype.push_back(haplotype);
							temp[t].count.push_back(depth);
						}
					}
					if (found == false)
					{
						temp.push_back(one_match());
						temp.back().cat_id = id;
						temp.back().count.push_back(depth);
						temp.back().haplotype.push_back(haplotype);
						temp.back().loc_id = loc_id;
					}
					
				}
				line_count++;
			}
			dad.close();
			for (size_t t = 0; t < temp.size(); t++)
			{
				if (temp[t].haplotype.size() == 2)//only keep those with 2 alleles
				{
					dad_off.dad_alleles.hap1.push_back(alleles());
					dad_off.dad_alleles.hap2.push_back(alleles());
					dad_off.off_alleles.hap1.push_back(alleles());
					dad_off.off_alleles.hap2.push_back(alleles());
					dad_off.dad_alleles.hap1[count].hap = temp[t].haplotype[0];
					dad_off.dad_alleles.hap2[count].hap = temp[t].haplotype[1];
					dad_off.dad_alleles.hap1[count].count = temp[t].count[0];
					dad_off.dad_alleles.hap2[count].count = temp[t].count[1];
					dad_off.dad_alleles.locus.push_back(temp[t].loc_id);
					dad_off.dad_alleles.catalog.push_back(temp[t].cat_id);
					dad_off.off_alleles.catalog.push_back(int());
					dad_off.off_alleles.locus.push_back(int());
					dad_off.mom_allele.push_back(string());
					count++;
				}
			}
			cout << "\nRead " << count << " loci from file " << dad_name << '\n';
			summary << '\t' << count;
			//now read the offspring file
			off.open(off_name);
			FileTest(off, off_name);
			line_count = count = 0;
			temp = vector<one_match>();
			temp.resize(0);
			while (universal_getline(off, line1))
			{
				if (line_count == 0)
				{
					stringstream ss1;
					ss1.str(line1);
					ss1 >> tmp1 >> tmp2 >> id >> samp_id >> loc_id >> haplotype >> depth >> log_like;
					temp.push_back(one_match());
					temp[0].cat_id = id;
					temp[0].count.push_back(depth);
					temp[0].haplotype.push_back(haplotype);
					temp[0].loc_id = loc_id;
				}
				else
				{
					stringstream ss1;
					ss1.str(line1);
					ss1 >> tmp1 >> tmp2 >> id >> samp_id >> loc_id >> haplotype >> depth >> log_like;
					//see if the offspring locus is in the dad's dataset<-match catalog ids
					for (size_t t = 0; t < dad_off.dad_alleles.catalog.size(); t++)
					{
						if (dad_off.dad_alleles.catalog[t] == id)
						{//need to save to temp to make sure there aren't more than 2 alleles
							bool found = false;
							for (size_t t = 0; t < temp.size(); t++)
							{
								if (temp[t].cat_id == id)
								{
									found = true;
									temp[t].haplotype.push_back(haplotype);
									temp[t].count.push_back(depth);
								}
							}
							if (found == false)
							{
								temp.push_back(one_match());
								temp.back().cat_id = id;
								temp.back().count.push_back(depth);
								temp.back().haplotype.push_back(haplotype);
								temp.back().loc_id = loc_id;
							}
						}
					}
				}
				line_count++;
			}//getline off
			off.close();
			count = 0;
			for (size_t t = 0; t < temp.size(); t++)
			{
				if (temp[t].haplotype.size() == 2)//only keep those with 1 or 2 alleles
				{
					for (size_t tt = 0; tt < dad_off.dad_alleles.catalog.size(); tt++)
					{
						if (dad_off.dad_alleles.catalog[tt] == temp[t].cat_id)
						{
							dad_off.off_alleles.hap1[tt].hap = temp[t].haplotype[0];
							dad_off.off_alleles.hap2[tt].hap = temp[t].haplotype[1];
							dad_off.off_alleles.hap1[tt].count = temp[t].count[0];
							dad_off.off_alleles.hap2[tt].count = temp[t].count[1];
							dad_off.off_alleles.locus[tt] = temp[t].loc_id;
							dad_off.off_alleles.catalog[tt] = temp[t].cat_id;
							count++;
						}
					}
				}
				if (temp[t].haplotype.size() == 1)//only keep those with 1 or 2 alleles
				{
					for (size_t tt = 0; tt < dad_off.dad_alleles.catalog.size(); tt++)
					{
						if (dad_off.dad_alleles.catalog[tt] == temp[t].cat_id)
						{
							dad_off.off_alleles.hap1[tt].hap = temp[t].haplotype[0];
							dad_off.off_alleles.hap2[tt].hap = temp[t].haplotype[0];
							dad_off.off_alleles.hap1[tt].count = temp[t].count[0];
							dad_off.off_alleles.hap2[tt].count = temp[t].count[0];
							dad_off.off_alleles.locus[tt] = temp[t].loc_id;
							dad_off.off_alleles.catalog[tt] = temp[t].cat_id;
							count++;
						}
					}
				}
			}
			cout << "\nFound " << count << " loci in offspring file " << off_name << " that matched parent 1 alleles.\n";
			summary << '\t' << count;
			//now infer female allele and write to file
			mom_name = out_path + "MOM" + dad_off.off_id.substr(3) + ".txt";
			mom.open(mom_name);
			output_name = out_path + dad_off.dad_id + "-" + dad_off.off_id + ".txt";
			output.open(output_name);
			cout << "\nWriting informative loci and inferred maternal alleles to file " << output_name << " and to " << mom_name << '\n';
			mom << "LocusID\tInferredAllele";
			output << "LocusID\tPar1Hap1\tPar1Cnt1\tPar1Hap2\tPar1Cnt2\tOffHap1\tOffCnt1\tOffHap2\tOffCnt2\tInferredPar2";
			for (size_t t = 0; t < dad_off.dad_alleles.catalog.size(); t++)
			{
				if (dad_off.off_alleles.catalog[t] == dad_off.dad_alleles.catalog[t]) //make sure they match
				{
					//compare dad's allele 1 to offspring alleles
					string mom_allele = "0";
					string dad1 = dad_off.dad_alleles.hap1[t].hap;
					string dad2 = dad_off.dad_alleles.hap2[t].hap;
					string kid1 = dad_off.off_alleles.hap1[t].hap;
					string kid2 = dad_off.off_alleles.hap2[t].hap;
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
						mom_allele = "0";
					if (dad1 == "0" || kid1 == "0")
						mom_allele = "0";
				
					dad_off.mom_allele[t] = mom_allele;
					output << '\n' << dad_off.dad_alleles.catalog[t] << '\t' << dad_off.dad_alleles.hap1[t].hap << '\t' << dad_off.dad_alleles.hap1[t].count
						<< '\t' << dad_off.dad_alleles.hap2[t].hap << '\t' << dad_off.dad_alleles.hap2[t].count << '\t' << dad_off.off_alleles.hap1[t].hap
						<< '\t' << dad_off.off_alleles.hap1[t].count << '\t' << dad_off.off_alleles.hap2[t].hap << '\t' << dad_off.off_alleles.hap2[t].count
						<< '\t' << dad_off.mom_allele[t];
					mom << '\n' << dad_off.dad_alleles.catalog[t] << '\t' << dad_off.mom_allele[t] << "\t0\t-\t0";
				}
			}
			output.close();
			mom.close();
			dad_name = off_name = mom_name = output_name = "";
			//update reference
			if (ii == 0)
			{//if it's the first one, need to create the consensus from the dad.
				for (size_t t = 0; t < dad_off.dad_alleles.locus.size(); t++)
				{
					consensus_loci.locus_id.push_back(dad_off.dad_alleles.locus[t]);
					consensus_loci.all_alleles.push_back(vector<string>());
					consensus_loci.all_alleles[t].push_back(dad_off.dad_alleles.hap1[t].hap);
					consensus_loci.all_alleles[t].push_back(dad_off.dad_alleles.hap2[t].hap);
					consensus_loci.num_in_dads.push_back(1);
					bool found1, found2;
					found1 = found2 = false;
					for (size_t tt = 0; tt < consensus_loci.all_alleles[t].size(); tt++)
					{
						if (dad_off.off_alleles.hap1[t].hap == consensus_loci.all_alleles[t][tt])
							found1 = true;
						if (dad_off.off_alleles.hap2[t].hap == consensus_loci.all_alleles[t][tt])
							found2 = true;
					}
					if (found1 == false)
						consensus_loci.all_alleles[t].push_back(dad_off.off_alleles.hap1[t].hap);
					if (found2 == false)
						consensus_loci.all_alleles[t].push_back(dad_off.off_alleles.hap2[t].hap);
					if (dad_off.mom_allele[t].size() > 0)
						consensus_loci.num_matched.push_back(1);
					else
						consensus_loci.num_matched.push_back(0);
				}
			}
			else
			{//after that we compare the individuals to the consensus.
				for (size_t t = 0; t < dad_off.dad_alleles.locus.size(); t++)
				{
					bool found_locus = false;
					for (size_t tt = 0; tt < consensus_loci.locus_id.size(); tt++)
					{
						if (consensus_loci.locus_id[tt] == dad_off.dad_alleles.locus[t])
						{//then it's the same locus
							found_locus = true;
							consensus_loci.num_in_dads[tt]++;
							if (dad_off.mom_allele[t].size() > 0)
								consensus_loci.num_matched[tt]++;
							bool found1, found2, found3, found4;
							found1 = found2 = found3 = found4 = false;
							for (size_t ttt = 0; ttt < consensus_loci.all_alleles[tt].size(); ttt++)
							{
								if (dad_off.dad_alleles.hap1[t].hap == consensus_loci.all_alleles[tt][ttt])
									found1 = true;
								if (dad_off.dad_alleles.hap2[t].hap == consensus_loci.all_alleles[tt][ttt])
									found2 = true;
								if (dad_off.off_alleles.hap1[t].hap == consensus_loci.all_alleles[tt][ttt])
									found3 = true;
								if (dad_off.off_alleles.hap2[t].hap == consensus_loci.all_alleles[tt][ttt])
									found4 = true;
							}
							if (found1 == false)
								consensus_loci.all_alleles[tt].push_back(dad_off.dad_alleles.hap1[t].hap);
							if (found2 == false)
								consensus_loci.all_alleles[tt].push_back(dad_off.dad_alleles.hap2[t].hap);
							if (found3 == false)
								consensus_loci.all_alleles[tt].push_back(dad_off.off_alleles.hap1[t].hap);
							if (found4 == false)
								consensus_loci.all_alleles[tt].push_back(dad_off.off_alleles.hap2[t].hap);
						}
					}
					if (found_locus == false)
					{ 
						int index = consensus_loci.locus_id.size();
						consensus_loci.locus_id.push_back(dad_off.dad_alleles.locus[t]);
						consensus_loci.all_alleles.push_back(vector<string>());
						consensus_loci.all_alleles[index].push_back(dad_off.dad_alleles.hap1[t].hap);
						consensus_loci.all_alleles[index].push_back(dad_off.dad_alleles.hap2[t].hap);
						consensus_loci.num_in_dads.push_back(1);
						bool found1, found2;
						found1 = found2 = false;
						for (size_t tt = 0; tt < consensus_loci.all_alleles[index].size(); tt++)
						{
							if (dad_off.off_alleles.hap1[t].hap == consensus_loci.all_alleles[index][tt])
								found1 = true;
							if (dad_off.off_alleles.hap2[t].hap == consensus_loci.all_alleles[index][tt])
								found2 = true;
						}
						if (found1 == false)
							consensus_loci.all_alleles[index].push_back(dad_off.off_alleles.hap1[t].hap);
						if (found2 == false)
							consensus_loci.all_alleles[index].push_back(dad_off.off_alleles.hap2[t].hap);
						if (dad_off.mom_allele[t].size() > 0)
							consensus_loci.num_matched.push_back(1);
						else
							consensus_loci.num_matched.push_back(0);
					}
				}//t
			}//ii > 0
			
			ii++;
		}//eof
	}//getline list
	paired_list.close();
	summary.close();
	cout << "\nWriting summary data per locus to file.\n";
	consensus_name = out_path + "loci_summary.txt";
	consensus.open(consensus_name);
	consensus << "LocusID\tNumAlleles\tAllelesCSV\tNumDadsWithLocus\tNumMomsInferred";
	for (size_t t = 0; t < consensus_loci.locus_id.size(); t++)
	{
		consensus << '\n' << consensus_loci.locus_id[t] << '\t' << consensus_loci.all_alleles.size() << '\t' << consensus_loci.all_alleles[t][0];
		for (size_t tt = 1; tt < consensus_loci.all_alleles[t].size(); tt++)
			consensus << "," << consensus_loci.all_alleles[t][tt];
		consensus << '\t' << consensus_loci.num_in_dads[t] << '\t' << consensus_loci.num_matched[t];
	}
	consensus.close();
	cout << "\nDone! Input integer to quit: \n";
	cin >> end;
	return 0;
}