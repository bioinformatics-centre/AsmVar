/************************
 * Author : Shujia Huang
 * Date   : 2014-08-01
 ************************/
#ifndef __VARUNIT_H__
#define __VARUNIT_H__

#include<iostream> 
#include<fstream>

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

    long score;
    double mismap; // mismatch probability

	// Reserve target. Can only use in recording translocations.
	// [Because Translocation have two target regions] 
    Region exp_target; 
    string exp_tarSeq;

public:

    VarUnit();
    VarUnit(const VarUnit &V);

    void ConvQryCoordinate(unsigned int qrySeqLen);
    void Swap();

	VarUnit ReAlign(Fa &target, Fa &query, AgeOption opt);

    void Clear(){ isClear = true; }
    bool Empty(){ return isClear; } // Do not output if isClear==true
	bool IsSuccessReAlign(){ return isSuccessReAlign; }

public:

    void OutStd(unsigned int tarSeqLen, unsigned int qrySeqLen, ofstream &O);
    void OutStd(unsigned int tarSeqLen, unsigned int exp_tarSeqLen, 
				unsigned int qrySeqLen, ofstream &O);
private:
	// Return the number of 'n' base in 'str'
	unsigned int NLength ( string &str ) {

    	unsigned int num(0);
    	for (size_t i(0); i < str.size(); ++i)
        	if (str[i] == 'N' || str[i] == 'n') ++num;
    	return num;
	}

private:
    bool isClear;
	bool isSuccessReAlign; // Check the variant could be re-align or not
};

/////////////////////////////////////////////////////////////////////////
/************************* Class AgeAlignment **************************/
/////////////////////////////////////////////////////////////////////////
class AgeAlignment {

public:

    AgeAlignment(): isInit_(false){}
    AgeAlignment(VarUnit &v, AgeOption opt){ Init(v, opt); }
	VarUnit vu() { return vu_; } 

    void Init(VarUnit &v, AgeOption opt);
    bool Align(string &tarFa, string &qryFa);

private:

    void ExtendVariant(unsigned long int, unsigned long int, int extandFlankSzie);
    bool IsHugeMemory(unsigned long int n, unsigned long int m);

private:

    VarUnit   vu_;
    AgeOption para_; // Parameters for AGE
    bool isInit_;
};

#endif



