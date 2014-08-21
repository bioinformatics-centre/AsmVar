/* 
 * Author : Shujia Huang
 * Date   : 2013-10-2
 *
 * This is the header of Variant.cpp
 *
 */
#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "utility.h"
#include "maf.h"
#include "VarUnit.h"
#include "MapReg.h"
#include "Region.h"
#include "Fa.h"
#include "vcf.h"

#include "AgeOption.h"

using namespace std;

// Rename 2014-06-19 14:38:02
class Variant : public MAF {

public: 

	Fa tarfa, qryfa;

private: 
	string sample;              // The name of sample
	vector< VarUnit > homoRef;  // The homozygous reference region, we should
						    	// filter the SV regions before output

	vector< VarUnit > nSeq;  // The reference or query is 'n' base or region
							 // it's juse the n base not in indel or any
							 // other SV region.

	vector< VarUnit > snp;           // Stored the SNP
	vector< VarUnit > intragap;      // Call the intra-scaffold-gap, just for the gaps which in alignment. Abort, 2014-02-27 19:28:42
	vector< VarUnit > insertion;     //
	vector< VarUnit > deletion;      //
	vector< VarUnit > inversion;     //
	vector< VarUnit > translocation; //
	vector< VarUnit > simulreg;// simultaneous gap regions
	vector< VarUnit > clipreg; // clip regions 
	vector< VarUnit > nomadic; // The query which totally can not map to target
	vector< VarUnit > nosolution;
	// Put all the VarUnit in this variant after 'AGE_Realign'
	map<string, vector<VarUnit> > allvariant; // refererce ID -> Variant

	// Make 'CallIversion' function to be private, because this function can 
	// only be call when all the alignments have been loaded in memerber value 
	// 'mapreg'.
	void CallReg        ( MapReg mapreg, string type, vector< VarUnit > & varReg );
	bool CallTranslocat ( MapReg left, MapReg middle, MapReg right );
	bool CallIversion   ( MapReg left, MapReg middle, MapReg right );
	bool CallSimultan   ( MapReg left, MapReg right );
	bool IsSameStrand   ( vector<MapReg> & mapreg );
	void FilterReg      ( map< string,vector<Region> >, map<string, size_t>, vector<VarUnit> & region );//just used in Filter ()
	void Output         ( vector< VarUnit > &, ofstream& O );
	void Assign2allvariant(vector<VarUnit> &v);
	void Unique(vector<VarUnit> &v); // Unique the variant in 'allvariant'
	void AGE_Realign(vector<VarUnit> &var); // AGE-Process

public :
	map< string, vector<MapReg> > mapreg;  // Mapping region. use for getting novel region
	map< string, vector<Region> > maptar;  // target's mapped regions
	map< string, vector<Region> > mapqry;  // query's  mapped regions
	map< string, long int > summary; // Used for record some summary information
	map< string, vector<Region> > VarTarRegs();// Return the target region of variants

public : // Can be called per-axt alignment. And will be called in main function
	void AssignSample  (string id) { sample = id; }
	void CallHomoRef   (); // call the homo region( the same as reference)
	void CallnSeq      (); // call the n region of target or query.
	void CallSNP       ();
	void CallInsertion ();
	void CallDeletion  ();
	void GetMapReg     ();
	void AGE_Realign   (); // AGE-Process

public : // Just can be call when all the axt alignments have been read!
	void CallSV        (); //It's SV not Indel, which cannot be called by a single alignment. simultaneous gaps,Inversion,translocation and no solution region 
	void CallClipReg   ();
	void CallNomadic   ();
	void Filter        ();  // Filter the variant regions which in nosolution
	void Output   ( string file ); // Output to the stdout
	void OutputSNP( string file ); // Output SNP
	void Summary  ( string file ); // Output Summary information
	void OutputGap( string file ); // Output the inter-gaps between different scaffold of the same target chromosome.

	void Output2VCF( string file ); // Output into a VCF format

	// Use for get the gap region, which actually would be the indel regions. can call deletion or insertion
	// Friend ship functions, but I fail to use friend function here, and I don't have enough time to figure out.
	vector< VarUnit > CallGap (Region& tar, string& tarSeq, Region& qry, string& qrySeq, char strand, long scroe, double mismap, string type);
	VarUnit CallGap ( MapReg left, MapReg right ); // call simultaneous gaps.
	long int Covlength ( vector<Region> mapreg );
};

///////////////////////////////////////////////////////////////////////////////
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
inline bool MySortByQryM ( MapReg  i, MapReg  j ) { 

	if ( i.query.id == j.query.id ) {
		return ( i.query.start  < j.query.start  ); 
	} else {
		return ( i.query.id < j.query.id );
	}
}
inline bool MySortByQryV ( VarUnit i, VarUnit j ) { 
	if ( i.query.id == j.query.id ) {
		return ( i.query.start  < j.query.start  ); 
	} else {
		return ( i.query.id < j.query.id );
	}
}
inline bool SortRegion( Region i, Region j){ 
	if ( i.id == j.id ) {
		return ( i.start  < j.start  ); 
	} else {
		return ( i.id < j.id );
	}
}

long int RegionMin( vector<Region> & region );
long int RegionMax( vector<Region> & region );
string ReverseAndComplementary( string & seq );

#endif

