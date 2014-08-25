/************************
 * Author : Shujia Huang
 * Date   : 2014-08-01
 ************************/
#ifndef __VARUNIT_H__
#define __VARUNIT_H__

#include<iostream> 
#include<fstream>
#include<stdlib.h>
#include<utility>
#include<vector>
#include<map>

#ifdef AGE_TIME
#include <sys/time.h>
#endif

//--- AGE includes ---
#include "AGEaligner.h"
#include "Sequence.h"
#include "Scorer.h"
#include "AgeOption.h" // alignement options

//--- Other include ---
#include "Region.h"
#include "Fa.h"        // VarUnit AgeAlign() function need Fa.h 
#include "utility.h"

using namespace std;

class AgeAlignment;

class VarUnit { 
// Well, it looks like the class 'Axt/MAF', but no! Totally different! 
// Use for record the varants.

public:
    Region target; // target or said reference
    Region query;  // the mapping fa

    string tarSeq; // Target sequence of variant
    string qrySeq; // before use this value we should consider the coversion 
				   // coordinate problem, to make sure we can get the correct 
				   // query sequece

    char   strand;
    string type;   // Variant type
	bool isHete;   // Heterozygousity or not

    long score;    // LAST Alignment score
    double mismap; // mismatch probability, could be used to calcute 
				   // the Mapping quality

	// This value will get after AGE process
	int homoRun;
	bool isSuccessAlign; // Can be re-align but the result may be not good
	bool isGoodReAlign;  // The result of realignment is good!
	pair<int, int> cipos, ciend; // Just for reference postion
    vector<pair<int,int> > identity; // Average, left, right

	// Reserve target. Can only use in recording translocations.
	// [Because Translocation have two target regions] 
    Region exp_target; 
    string exp_tarSeq;

public:

    VarUnit();
    VarUnit(const VarUnit &V);

    void ConvQryCoordinate(long int qrySeqLen);
    void Swap();

	// Here is the only function that will use 'AgeAlignment' class
	vector<VarUnit> ReAlignAndReCallVar(string &target, string &query, 
										AgeOption opt);

    void Clear(){ tarSeq.clear(); qrySeq.clear(); isClear = true; }
    bool Empty(){ return isClear; } // Do not output if isClear==true

public:

	void OutErr();
    void OutStd(long int tarSeqLen, long int qrySeqLen, ofstream &O);
    void OutStd(long int tarSeqLen, long int exp_tarSeqLen, 
				long int qrySeqLen, ofstream &O);
private:
	// Return the number of 'n' base in 'str'
	long int NLength ( string &str ) {

    	long int num(0);
    	for (size_t i(0); i < str.size(); ++i)
        	if (str[i] == 'N' || str[i] == 'n') ++num;
    	return num;
	}

private:
    bool isClear;
};

vector<VarUnit> MergeVarUnit(vector<VarUnit> &VarUnitVector, int delta);

/////////////////////////////////////////////////////////////////////////
/************************* Class AgeAlignment **************************/
/////////////////////////////////////////////////////////////////////////
// Each single Indel or SV should AGE re-alignement. That's why I keep 
// 'AgeAlignment' class together with 'VarUnit' class!
// 'AgeAlignment' will just present in the memeber function:'ReAlign...()'
// of 'VarUnit' class! 
class AgeAlignment {

public:

    AgeAlignment(): isInit_(false), isalign_(false), isgoodAlign_(false) {}
    AgeAlignment(VarUnit &v, AgeOption opt) { Init(v, opt); }
	VarUnit vu() { return vu_; } 

    void Init(VarUnit &v, AgeOption opt);
    bool Align(string &tarFa, string &qryFa);

	// Caution: The position order is big->small in query if 
	// mapping to the '-' strand! We should take care about this!
	AlignResult AligResult() { return alignResult_; }	
	vector<VarUnit> VarReCall();

	bool isgoodAlign() { return isgoodAlign_; }
	int HomoRun() { return alignResult_._homo_run_atbp1; }
	pair<int, int> cipos() { return alignResult_._ci_start1; }
	pair<int, int> ciend() { return alignResult_._ci_end1;   }

private:

	// Variant in excise region
	VarUnit CallVarInExcise(pair<MapData, MapData> &m1, 
							pair<MapData, MapData> &m2, char strand);
	vector<VarUnit> CallVarInFlank(pair<MapData, MapData> &m, string &mapInfo,
								   char strand);
    void ExtendVU(long int, long int, int extandFlankSzie);
    bool IsHugeMemory(long int n, long int m);

private:

    VarUnit   vu_;   // this would aways be a copy of a VarUnit.
    AgeOption para_; // Parameters for AGE
	AlignResult alignResult_;
    bool isInit_;
	bool isalign_;
	bool isgoodAlign_;
};


#endif



