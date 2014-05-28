/*
 * Author : Shujia Huang
 * Date   : 2013-11-23 00:38:59
 *
 * */
#ifndef STATS_H
#define STATS_H

#include <iostream>
#include <vector>
#include <map>

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "tabix/tabix.hpp"
#include "utility.h"    // My own header file
#include "Coverage.h"

using namespace std;
using namespace BamTools;

const string ERROR_PREFIX = "[Main stats] ";

const short INNER      = 1;
const short OUTTER     = 2;
const short SAME       = 3;
const short DIFF_CHROM = 4;
const short UNPAIRED   = 5;

// I should build Stats to be a class -_-!
struct Stats {

	map <string, pair<int,int> > insert;
	vector<int> refGC;
    vector<float> gc2cov;

	Coverage properPair;     // This is not as the proper pair as bwa. We consider the aligment score at the same time
	// UnproperPair : Include lowAlignScore,badInsertSizePair,WrongOrientation,
	Coverage pLowAlignScore; // This is proper pair by bwa but with low aligment score.
	Coverage crossSignal;    // The reads contain 'I' signal or 'D' signal in cigar string in bam file
	Coverage badInserSizePair;
	Coverage wrongOrientation;
	Coverage single;         // Single-end reads
	Coverage clipAlign;

	// in sorted BAM file, the reads are in order, which means the fragments, where
	// we include the reads as part of the fragments, are ordered by start position.
	// But they are not in order if we don't include the reads, since reads could
	// be of different lengths or clipped. Use a multiset so that fragments are automatically
	// kept in order
};

short GetPairOrientation(BamAlignment& al, bool reverse=false);

//record the 'I' and 'D' signal by reference coordinates in the cigar string!
map< char, vector<pair<int32_t,int32_t> > > IndelSignalRefPos( int32_t pos, vector<CigarOp> CigarData );

// Loads the contents of file into vector.  The file is expected to be
// of the form:
// GC      coverage
// (tab separated), where GC is all integers in [0,100] in ascending numerical order.
// Used to be the expected coverage at each position
void LoadGC2Cov(string& filename, vector<float>& v_in);

// Loads the GC from bgzipped tabixed file, for just the given reference sequence.
// Fills the vector with GC content at each position of the sequence (v_in[n] = GC at position n (zero based))
void LoadGC(Tabix& ti, string& refname, vector<int>& v_in);

// fills each list with (start, end) positions of gaps in the input file, of the form:
// sequence_name<TAB>start<TAB>end
// and this file is expected to be bgzipped and tabixed
// Map key = sequence_name
void LoadGaps(string fname, map<string, list<pair<long,long> > >& gaps);

// The start coordinate of coverage should be the same in stats if the position isn't 0 
bool CheckStatsAndFitLength ( Stats& stats );

// Out put all the covered stats of each position left to 'pos', not include 'pos'
void OutputStats (string refname, long pos, Stats& stats, map<string, list<pair<long, long> > >& gaps);

#endif

