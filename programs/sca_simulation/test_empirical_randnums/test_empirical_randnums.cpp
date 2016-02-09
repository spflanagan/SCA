#pragma once
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

// The following code has to do with the random number
// generator. This is the Mersenne Twister prng.
// The reference is:
// M. Matsumoto and T. Nishimura,
// Mersenne Twister: A 623-Dimensionally Equidistributed
// Uniform Pseudo-Random Number Generator,
// ACM Transactions on Modeling and Computer Simulation,
// Vol. 8, No. 1, January 1998, pp. 3-30
// Copyright corresponding to Mersenne Twister Code:
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. The names of its contributors may not be used to endorse or promote
// products derived from this software without specific prior written
// permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)
#define My_PI 3.14159265358979323846

static unsigned long mt[N]; /* the array for the state vector */
static int mti = N + 1; /* mti==N+1 means mt[N] is not initialized */
double Two2the36 = 4294967296.0;
static double rnZ = 0;

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long int seed)
{
	/* setting initial seeds to mt[N] using the generator Line 25 of Table 1 in	[KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102] */
	mt[0] = seed & 0xffffffff;
	for (mti = 1; mti<N; mti++)
		mt[mti] = (69069 * mt[mti - 1]) & 0xffffffff;
}

inline double genrand()
{
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	/* mag01[x] = x * MATRIX_A for x=0,1 */
	if (mti >= N) { /* generate N words at one time */
		int kk;
		if (mti == N + 1) /* if sgenrand() has not been called, */
			sgenrand(4357); /* a default initial seed is used */
		for (kk = 0; kk<N - M; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1];
		}

		for (; kk<N - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U(y);
	y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L(y);

	return ((double)y / Two2the36); /* reals */
	//return y; /* for integer generation */
	// This should give a range of [0,1); for [0,1]
	// use Two2the36-1.
}

inline int randnum(int Max) // returns a value between 0 and Max-1.
{
	double rnms;
	double dbrnum;
	int irnum;
	rnms = genrand();
	dbrnum = floor(rnms * Max);
	irnum = dbrnum;
	return irnum;
}

//this function was written by Sarah P. Flanagan
//17 December 2015
//It reads in a file containing known allele frequency distributions
//and returns an integer for which allele index value to use
using namespace std;
int d_empirical_afs(int num_alleles, string afs_file_name)
{
	//read in the known distribution
	ifstream afs_file;
	string line;
	vector<double> intervals;
	vector<double> densities;
	vector<int> count_intervals;
	double t, tt, last_interval, last_prob, sum;
	int j, jj, jjj, int_index, count, return_allele;
	bool found;

	afs_file.open(afs_file_name);
	while (getline(afs_file, line))
	{
		stringstream temp;
		temp.str(line);
		temp >> t >> tt;
		intervals.push_back(t);
		densities.push_back(tt);
	}
	//need to divide up the known distribution for the correct number of alleles
	vector<double> allele_probs;
	int num_allele_intervals = intervals.size() / num_alleles;
	if (num_allele_intervals*num_alleles < intervals.size())
		num_allele_intervals++;
	for (j = 0; j < num_alleles; j++)
		allele_probs.push_back(0);
	int_index = count = 0;
	count_intervals.push_back(0);
	for (j = 0; j < intervals.size(); j++)
	{
		if (count >= num_allele_intervals)
		{
			count_intervals.push_back(0);
			int_index++;
			count = 0;
		}
		allele_probs[int_index] = allele_probs[int_index] + densities[j];
		count_intervals[int_index]++;
		count++;
	}
	sum = 0;
	for (j = 0; j < num_alleles; j++)
	{
		sum = sum + allele_probs[j];
		allele_probs[j] = allele_probs[j] / count_intervals[j];
	}
	for (j = 0; j < num_alleles; j++)
		allele_probs[j] = allele_probs[j] / sum;

	//now need to draw a number from that distribution and know which allele it should be.
	found = false;
	while (!found)
	{
		count = randnum(num_alleles);
		t = genrand();
		if (count*t < allele_probs[count])
		{
			return_allele = count;
			found = true;
		}
	}
	return return_allele;
}

int main()
{
	int i, ii, num_alleles;
	string out_name, in_name;
	ofstream out;

	num_alleles = 6;
	in_name = "E://ubuntushare//SCA//sca_simulation//allelefreqs.txt";
	out_name = "test_distribution.txt";
	out.open(out_name);
	for (i = 0; i < 10000; i++)
	{
		ii = d_empirical_afs(num_alleles, in_name);
		out << '\n' << ii;
	}
	out.close();
	cout << "\nDone! Input integer to quit.\n";
	cin >> i;
	return 0;
}