/* 
 * Author : Shujia Huang
 * Date   : 2013-10-2
 *
 * This is the header of Variant.cpp
 *
 */
#ifndef VARIANT_H
#define VARIANT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "maf.h"
#include "Region.h"
#include "utility.h"

using namespace std;

unsigned int NLength ( string & str ); // Used to calculate the 'n' length in class 'VarUnit' and 'Variant'

class VarUnit { // Well, it looks like the class 'Axt/MAF', but no! Totally different! Use for record the varants.

private:
	bool isClear;

public:
	Region target; // target or said reference
	Region query;  // the mapping one

	string tarSeq;
	string qrySeq; // if include this part we should consider the coversion coordinate problem, so that make me to get the right seq
	char   strand;
	string type;	

	long score;
    double mismap; // mismatch probability

	Region exp_target; // Reserve target. Can only use in recording translocations. [Because Translocation have two target regions]
	string exp_tarSeq;

public:
	VarUnit () { 
		target.id = "-"; query.id = "-"; tarSeq = "-"; qrySeq = "-"; strand = '.', type = "." ; isClear = false; 
		score     = 0;   mismap   = 1.0;
	}
	VarUnit ( const VarUnit & V ) { 
		target = V.target; query = V.query; tarSeq = V.tarSeq; qrySeq = V.qrySeq; strand = V.strand; 
		type   = V.type; isClear = V.isClear; score= V.score;  mismap = V.mismap;
		exp_target = V.exp_target;
	}

	void ConvQryCoordinate ( unsigned int qrySeqLen ) { 
	// This funtion just conversion the coordinate of Axt/MAF format creat by 'lastz'/'last ', which mapped to the '-' strand
		if ( strand != '-' ) return;
		unsigned int itemp = query.start;
        query.start = qrySeqLen - query.end + 1;
        query.end   = qrySeqLen - itemp + 1;
	}

	void Swap () { // Swap the target and query region. Ignore 'exp_target' 
		Region tmp = target; target = query;  query  = tmp;
		string str = tarSeq; tarSeq = qrySeq; qrySeq = str;
	}
	void Clear() { isClear = true; }
	bool Empty() { return isClear; } // Do not output if isClear==true

public:
	void OutStd ( unsigned int tarSeqLen, unsigned int qrySeqLen, ofstream& O ) { // Output the axt alignment to STDERR

		if ( tarSeq.empty() || qrySeq.empty() ) { cerr << "tarSeq.empty() || qrySeq.empty()" << endl; exit(1); }

		unsigned int qnl = NLength ( qrySeq );
		unsigned int tnl = NLength ( tarSeq );
		O   << target.id << "\t" << target.start << "\t" << target.end << "\t" << target.end - target.start + 1    << "\t"
			<< double(tnl)/tarSeq.length()       << "\t" << tarSeqLen  << "\t" << query.id  << "\t" << query.start << "\t" 
			<< query.end << "\t" << query.end  - query.start  + 1      << "\t" << double(qnl)/qrySeq.length()      << "\t"
			<< qrySeqLen << "\t" << strand       << "\t" << score      << "\t" << mismap    << "\t" << type        << endl;
	}
	void OutStd ( unsigned int tarSeqLen, unsigned int exp_tarSeqLen, unsigned int qrySeqLen, ofstream& O ) {
		if ( exp_target.isEmpty() ) cerr << "[ERROR]exp_target is empty!\n";

		OutStd ( tarSeqLen, qrySeqLen, O );

		if ( exp_tarSeq.empty() ) { cerr << "exp_tarSeq.empty() " << endl; exit(1); }

        unsigned int qnl = NLength ( qrySeq );
        unsigned int tnl = NLength ( exp_tarSeq );
        O   << exp_target.id << "\t" << exp_target.start << "\t" << exp_target.end << "\t" << exp_target.end - exp_target.start + 1 << "\t"
            << double(tnl)/exp_tarSeq.length()           << "\t" << exp_tarSeqLen  << "\t" << query.id  << "\t"   << query.start    << "\t"
            << query.end << "\t" << query.end  - query.start  + 1          << "\t" << double(qnl)/qrySeq.length() << "\t"
            << qrySeqLen << "\t" << strand <<"\t"<< score<< "\t" << mismap << "\t" << type + "-E" << endl;
	}
};

class MapReg { // Mapping regions

public:
	Region target; // target or said reference
	Region query;  // the mapping one
	char  strand;  // Mapping strand

    long score;
    double mismap; // mismatch probability

public: 
	void OutErrReg() {
		cerr << "# " << score    << "\t" << mismap << "\t" << target.id << "\t" << target.start << "\t" << target.end 
			 << "\t" << query.id << "\t" << query.start  << "\t" << query.end  << "\t" << strand << endl;
	}
};

// Rename 2014-06-19 14:38:02
class Variant : public MAF {

private: 
	string sample;                    // The name of sample
	vector< VarUnit > snp;           // Stored the SNP
	vector< VarUnit > intragap;      // Call the intra-scaffold-gap, just for the gaps which in alignment. Abort, 2014-02-27 19:28:42
	vector< VarUnit > insertion;     //
	vector< VarUnit > deletion;      //
	vector< VarUnit > inversion;     //
	vector< VarUnit > translocation; //
	vector< VarUnit > simulreg;      // simultaneous gap regions
	vector< VarUnit > clipreg;       // clip regions 
	vector< VarUnit > nomadic;       // The query which totally can not map to target
	vector< VarUnit > nosolution;

	// Make 'CallIversion' function to be private, because this function can only be call when all
	// the alignments have been loaded in memerber value 'mapreg'.
	void CallReg        ( MapReg mapreg, string type, vector< VarUnit > & varReg );
	bool CallTranslocat ( MapReg left, MapReg middle, MapReg right );
	bool CallIversion   ( MapReg left, MapReg middle, MapReg right );
	bool CallSimultan   ( MapReg left, MapReg right );
	//void CallNoSolution ( MapReg mapreg );
	bool IsSameStrand   ( vector<MapReg> & mapreg );
	void FilterReg      ( map< string,vector<Region> >, map<string, size_t>, vector<VarUnit> & region );//just used in Filter ()
	void Output         ( vector< VarUnit > &, ofstream& O );

public :
	map< string, vector<MapReg> > mapreg;  // Mapping region. use for getting novel region
	map< string, vector<Region> > maptar;  // target's mapped regions
	map< string, vector<Region> > mapqry;  // query's  mapped regions
	map< string, unsigned int   > summary; // Used for record some summary information
	map< string, vector<Region> > VarTarRegs();// Return the target region of variants

public : // Can be called per-axt alignment. And will be called in main function
	void AssignSample  (string id) { sample = id; }
	void CallSNP       ();
	void CallInsertion ();
	void CallDeletion  ();
	void GetMapReg     ();

public : // Just can be call when all the axt alignments have been read!
	void CallSV        (); //It's SV not Indel, which cannot be called by a single alignment. simultaneous gaps,Inversion,translocation and no solution region 
	void CallClipReg   ();
	void CallNomadic   ();
	void Filter        ();  // Filter the variant regions which in nosolution
	void Output   ( string file ); // Output to the stdout
	void OutputSNP( string file ); // Output SNP
	void Summary  ( string file ); // Output Summary information
	void OutputGap( string file ); // Output the intra-gap in the scaffold and the inter-gaps between different scaffold of the same target chromosome.

	// Use for get the gap region, which actually would be the indel regions. can call deletion or insertion
	// Friend ship functions, but I fail to use friend function here, and I don't have enough time to figure out.
	vector< VarUnit > CallGap (Region& tar, string& tarSeq, Region& qry, string& qrySeq, char strand, long scroe, double mismap, string type);
	VarUnit CallGap ( MapReg left, MapReg right ); // call simultaneous gaps.
	unsigned int Covlength ( vector<Region> mapreg );
};
//inline bool MySortByTarM ( MapReg  i, MapReg  j );
//inline bool MySortByTarV ( VarUnit i, VarUnit j );
//inline bool MySortByQryM ( MapReg  i, MapReg  j );
//inline bool MySortByQryV ( VarUnit i, VarUnit j );
//inline bool SortRegion   ( Region  i, Region  j );
inline bool MySortByTarM ( MapReg  i, MapReg  j ) {
    if ( i.target.id == j.target.id ) {
        return ( i.target.start < j.target.start );
    } else {
        return ( i.target.id <  j.target.id );
    }
}
inline bool MySortByTarV ( VarUnit i, VarUnit j ) {
    if ( i.target.id == j.target.id ) {
        return ( i.target.start < j.target.start );
    } else {
        return ( i.target.id < j.target.id );
    }
}
inline bool MySortByQryM ( MapReg  i, MapReg  j ) { return ( i.query.start  < j.query.start  ); }
inline bool MySortByQryV ( VarUnit i, VarUnit j ) { return ( i.query.start  < j.query.start  ); }
inline bool SortRegion   ( Region  i, Region  j ) { return ( i.start < j.start  ); }

unsigned int RegionMin   ( vector<Region> & region );
unsigned int RegionMax   ( vector<Region> & region );
string ReverseAndComplementary ( string & seq );

#endif

