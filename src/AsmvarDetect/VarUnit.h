/************************
 * Author : Shujia Huang
 * Date   : 2014-08-01
 ************************/
#ifndef __VARUNIT_H__
#define __VARUNIT_H__

#include<iostream> 
#include<fstream>

#include "Region.h"

class VarUnit { // Well, it looks like the class 'Axt/MAF', but no! Totally different! Use for record the varants.

public:
    Region target; // target or said reference
    Region query;  // the mapping fa

    string tarSeq;
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

    VarUnit () {
        target.id = "-"; 
		query.id  = "-"; 
		tarSeq    = "."; 
		qrySeq    = "."; 
		strand    = '.'; 
		type      = "."; 
        score     = 0;   
		mismap    = 1.0;
		isClear   = false;
    }

    VarUnit ( const VarUnit& V ) {

        target = V.target; query   = V.query; tarSeq = V.tarSeq; 
		qrySeq = V.qrySeq; strand  = V.strand;
        type   = V.type  ; isClear = V.isClear; 
		score  = V.score ;  mismap = V.mismap;
        exp_target = V.exp_target;
		exp_tarSeq = V.exp_tarSeq;
    }

    void ConvQryCoordinate ( unsigned int qrySeqLen ) {
    // This funtion just conversion the coordinate of Axt/MAF format creat 
    // by 'lastz'/'last ', which mapped to the '-' strand
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

    void OutStd ( unsigned int tarSeqLen, unsigned int qrySeqLen, ofstream &O ) { // Output the axt alignment to STDERR

        if ( tarSeq.empty() || qrySeq.empty() ) { std::cerr << "tarSeq.empty() || qrySeq.empty()" << endl; exit(1); }

        unsigned int qnl = NLength ( qrySeq );
        unsigned int tnl = NLength ( tarSeq );
        O << target.id << "\t" << target.start << "\t" << target.end << "\t" << target.end - target.start + 1    << "\t"
          << double(tnl)/tarSeq.length()       << "\t" << tarSeqLen  << "\t" << query.id  << "\t" << query.start << "\t"
          << query.end << "\t" << query.end  - query.start  + 1      << "\t" << double(qnl)/qrySeq.length()      << "\t"
          << qrySeqLen << "\t" << strand       << "\t" << score      << "\t" << mismap    << "\t" << type        << endl;
    }

    void OutStd ( unsigned int tarSeqLen, unsigned int exp_tarSeqLen, unsigned int qrySeqLen, ofstream &O ) {
        if ( exp_target.isEmpty() ) cerr << "[ERROR]exp_target is empty!\n";

        OutStd ( tarSeqLen, qrySeqLen, O );

        if ( exp_tarSeq.empty() ) { cerr << "exp_tarSeq.empty() " << endl; exit(1); }

        unsigned int qnl = NLength ( qrySeq );
        unsigned int tnl = NLength ( exp_tarSeq );
        O << exp_target.id << "\t" << exp_target.start << "\t" << exp_target.end << "\t" << exp_target.end - exp_target.start + 1 << "\t"
          << double(tnl)/exp_tarSeq.length()           << "\t" << exp_tarSeqLen  << "\t" << query.id  << "\t"   << query.start    << "\t"
          << query.end << "\t" << query.end  - query.start  + 1          << "\t" << double(qnl)/qrySeq.length() << "\t"
          << qrySeqLen << "\t" << strand <<"\t"<< score<< "\t" << mismap << "\t" << type + "-E" << endl;
    }

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

};

#endif



