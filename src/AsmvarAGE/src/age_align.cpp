/****************************************************************************
 *
 *     AGE -- Alignment with Gap Excision
 *     Copyright (C)  Alexej Abyzov & Shujia Huang
 *                
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Creative Commons license
 * (Attribution-NonCommerical).
 * See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
 *
 * Author: Alexej Abyzov
 * Rewrite by : Shujia Huang 2014-4-8 16:26:53
 *
 */

//--- C/C++ includes ---
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <vector>
#include <string>

#ifdef AGE_TIME
#include <sys/time.h>
#endif

using namespace std;

//--- AGE includes ---
#include "AGEaligner.hh"
#include "Sequence.hh"
#include "Scorer.hh"
#include "gzstream.h"

//--- My own lib ---
#include "Region.h"
#include "Fa.h"
#include "Join.h"

class CmdLineOptions {

public :
	int gapOpen;
	int gapExtend;
	int match;
	int mismatch;
	//string coor1;
	//string coor2;
	string varInfile, tarFaInfile, qryFaInfile, selectTarId;
	bool indel, inv, invl, invr, tdup;
	bool both, revcom1, revcom2;
public :
	CmdLineOptions() : match(1), mismatch(-2), gapOpen(-2), gapExtend(-1), 
					   both(false) , revcom1(false), revcom2(false),
					   indel(false), inv(false), invl(false), invr(false),
					   tdup(false)
	{ selectTarId = "ALL"; }
};

struct VarReg {
    Region target; // target or said reference
    Region query;  // the mapping one
	string info;
};

void Usage ( const char* prog, CmdLineOptions& ops ) {

	cerr << "\nAGE -- Alignment with Gap Excision\n\n";
	cerr << "Author : Alexej Abyzov\nRewrite by : Shujia Huang on 2014-4-8 16:26:53\n\n";
	cerr << "Usage : " << prog << " [Options] [VariantRegionInfile] > Output.age 2> alig.log \n\n";
	cerr << "      Options  : \n";
	cerr << "          -t, --target    [str]  target sequence, fa format. require!\n";
	cerr << "          -q, --query     [str]  Query  sequence, fa format. require!\n";
	cerr << "          -s, --selectId  [str]  Select the target id. [" << ops.selectTarId << "]\n";
	cerr << "          -m, --match     [int]  match score.          [" << ops.match       << "]\n";
	cerr << "          -M, --mismatch  [int]  mismatch score.       [" << ops.mismatch    << "]\n";
	cerr << "          -g, --go        [int]  Gap open score.       [" << ops.gapOpen     << "]\n";
	cerr << "          -e, --ge        [int]  Gap Extend score.     [" << ops.gapExtend   << "]\n";
	cerr << "          -b, --both             Both reverse.         [" << ops.both        << "]\n";
	cerr << "          -a, --revcom1          revcom1.              [" << ops.revcom1     << "]\n";
	cerr << "          -c, --revcom2          revcom2.              [" << ops.revcom2     << "]\n";
	cerr << "          -i, --indel            Indel.                [" << ops.indel       << "]\n";
	cerr << "          -v, --inv              Inv.                  [" << ops.inv         << "]\n";
	cerr << "          -l, --invl             Invl.                 [" << ops.invl        << "]\n";
	cerr << "          -d, --tdup             tdup.                 [" << ops.tdup        << "]\n";
	cerr << "          -r, --invr             Invr.                 [" << ops.invr        << "]\n";
	cerr << "          -h                    Output this help information.\n\n";
	
	exit(1);
}

void ParseOptions ( int argc, char** argv, CmdLineOptions& ops ) {

	struct option longopts[] = {
		{"go"      , 1, NULL, 'g'},
		{"ge"      , 1, NULL, 'e'},
		{"match"   , 1, NULL, 'm'},
		{"mismatch", 1, NULL, 'M'},
		{"indel"   , 0, NULL, 'i'},
		{"inv"     , 0, NULL, 'v'},
		{"invl"    , 0, NULL, 'l'},
		{"invr"    , 0, NULL, 'r'},
		{"tdup"    , 0, NULL, 'd'},
		{"revcom1" , 0, NULL, 'a'},
		{"revcom2" , 0, NULL, 'c'},
		{"both"    , 0, NULL, 'b'},
		{"target"  , 1, NULL, 't'}, // Target targeterence in axt result. [fa format]
		{"query"   , 1, NULL, 'q'}, // Query  targeterence in axt result. [fa format]
		{"selectId", 1, NULL, 's'}, // Select the target id
		{0,0,0,0}
	};

	char c;
	while ( (c = getopt_long(argc, argv, "g:e:m:M:t:q:s:abcdhilrv", longopts, NULL) ) != -1 ) {
		switch ( c ) {
			case 'g' : ops.gapOpen     = atoi( optarg ); break;
			case 'e' : ops.gapExtend   = atoi( optarg ); break;
			case 'm' : ops.match       = atoi( optarg ); break;
			case 'M' : ops.mismatch    = atoi( optarg ); break;
			case 't' : ops.tarFaInfile = optarg        ; break;
			case 'q' : ops.qryFaInfile = optarg        ; break;
			case 's' : ops.selectTarId = optarg        ; break;
			case 'a' : ops.revcom1     = true          ; break;
			case 'c' : ops.revcom2     = true          ; break;
			case 'b' : ops.both        = true          ; break;
			case 'i' : ops.indel       = true          ; break;
			case 'v' : ops.inv         = true          ; break;
			case 'l' : ops.invl        = true          ; break;
			case 'r' : ops.invr        = true          ; break;
			case 'd' : ops.tdup        = true          ; break;
			case 'h' : Usage( argv[0], ops );

			default :
				cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl; exit(1);
		}
	}
	if (argc > 1) ops.varInfile = argv[argc-1];
	if (ops.varInfile.empty() )  { cerr << "\n[ERROR] Variant in file is empty!\n"; Usage(argv[0], ops); }
	if (ops.tarFaInfile.empty()) { cerr << "\n[ERROR] Target Fa is required!   \n"; Usage(argv[0], ops); }
	if (ops.qryFaInfile.empty()) { cerr << "\n[ERROR] Query  Fa is required!   \n"; Usage(argv[0], ops); }
	cerr << "\n# [INFO] The parameter of this program: " << join (argv, argc) << "\n" << endl;
    cerr << "    [INFO] Options  : \n";
    cerr << "          -t, --target    [" << ops.tarFaInfile << "]\n";
    cerr << "          -q, --query     [" << ops.qryFaInfile << "]\n";
    cerr << "          -s, --selectId  [" << ops.selectTarId << "]\n";
    cerr << "          -m, --match     [" << ops.match       << "]\n";
    cerr << "          -M, --mismatch  [" << ops.mismatch    << "]\n";
    cerr << "          -g, --go        [" << ops.gapOpen     << "]\n";
    cerr << "          -e, --ge        [" << ops.gapExtend   << "]\n";
    cerr << "          -b, --both      [" << ops.both        << "]\n";
    cerr << "          -a, --revcom1   [" << ops.revcom1     << "]\n";
    cerr << "          -c, --revcom2   [" << ops.revcom2     << "]\n";
    cerr << "          -i, --indel     [" << ops.indel       << "]\n";
    cerr << "          -v, --inv       [" << ops.inv         << "]\n";
    cerr << "          -l, --invl      [" << ops.invl        << "]\n";
    cerr << "          -d, --tdup      [" << ops.tdup        << "]\n";
    cerr << "          -r, --invr      [" << ops.invr        << "]\n\n";
	return;
}


vector <VarReg> LoadVariantRegion ( const char* file ) {

	vector< VarReg > varReg;
	igzstream I ( file );
    if ( !I ) { cerr << "# [ERROR] Cannot open file : " << file << endl; exit(1); }
    string tmp;
	VarReg reg;
    while ( 1 ) {
	// Format: TargetId  TargetStart  TargetEnd  QueryId  QueryStart QueryEnd  Info
        I >> reg.target.id;
        if ( I.eof() ) break;
		if (reg.target.id[0] == '#') { getline( I, tmp, '\n'); continue; }

		I >> reg.target.start >> reg.target.end >> reg.query.id >> reg.query.start >> reg.query.end >> reg.info;
        getline( I, tmp, '\n');

		varReg.push_back( reg );
    }
    I.close();

    return varReg;
}

////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[]) {

	CmdLineOptions ops;
    ParseOptions(argc, argv, ops);

	int flag = 0;
	if (ops.indel) flag |= AGEaligner::INDEL_FLAG;
	if (ops.inv  ) flag |= AGEaligner::INVERSION_FLAG;
	if (ops.invl ) flag |= AGEaligner::INVL_FLAG;
	if (ops.invr ) flag |= AGEaligner::INVR_FLAG;
	if (ops.tdup ) flag |= AGEaligner::TDUPLICATION_FLAG;
	if (flag == 0) flag  = AGEaligner::INDEL_FLAG; // Default

	int n_modes = 0;
    if (flag & AGEaligner::INDEL_FLAG)        n_modes++;
    if (flag & AGEaligner::TDUPLICATION_FLAG) n_modes++;
    if (flag & AGEaligner::INVR_FLAG)         n_modes++;
    if (flag & AGEaligner::INVL_FLAG)         n_modes++;
    if (flag & AGEaligner::INVERSION_FLAG)    n_modes++;

	bool error(false);
	if (n_modes != 1) {
		cerr << "# [Error] In mode specification. ";
		if ( n_modes == 0 ) cerr << "No mode is specified.\n";
		if ( n_modes  > 1 ) cerr << "More than one mode is specified.\n";
		error = true;
	}
	if ( error ) exit(1);

	Fa tarSeq, qrySeq;
	tarSeq.Load( ops.tarFaInfile );
	qrySeq.Load( ops.qryFaInfile );
	cerr << "# [INFO] Target fa and Query fa are loaded successfully." << local_time();
	
	vector<VarReg> varReg = LoadVariantRegion ( ops.varInfile.c_str() );
	cerr << "# [INFO] Variant regions are loaded successfully."        << local_time();
	if  (varReg.size() == 0) cerr << "[WARNING] The variant list is empty!\n";
	for ( size_t i(0); i < varReg.size(); ++i ) {

		if ((toupper(ops.selectTarId) != "ALL") && (toupper(ops.selectTarId) != toupper(varReg[i].target.id))) 
			continue;

		if ( !tarSeq.fa.count(varReg[i].target.id) ) { 
			cerr << "# [ERROR] " << varReg[i].target.id << " is not in " << ops.tarFaInfile << "\n";
			exit(1);
		}
		if ( !qrySeq.fa.count(varReg[i].query.id) ) { 
            cerr << "# [ERROR] " << varReg[i].query.id << " is not in " << ops.qryFaInfile << "\n";
            exit(1);
        }
		Sequence* old_seq1 = new Sequence( tarSeq.fa[varReg[i].target.id], varReg[i].target.id, 1, false );
		Sequence* old_seq2 = new Sequence( qrySeq.fa[varReg[i].query.id] , varReg[i].query.id , 1, false );
		
		cout << "# " << varReg[i].info << "\n";
		for (Sequence* s1 = old_seq1; s1; s1 = s1->next() ) {
			for (Sequence* s2 = old_seq2; s2; s2 = s2->next()) {
				Sequence* seq1 = s1->substr(varReg[i].target.start, varReg[i].target.end);
				Sequence* seq2 = s2->substr(varReg[i].query.start , varReg[i].query.end );

				if (ops.revcom1) seq1->revcom();
				if (ops.revcom2) seq2->revcom();
				// Sequence seq1("ACTGGTGTCAACTG"),seq2("ACTGACTG");
				// Sequence seq1("GGACTGGACTG"),   seq2("CACTGGACTG");
				// Sequence seq1("CCCCAGCTTT"),    seq2("AGCCTTT");
				// Sequence seq1("AGCCTTT"),       seq2("AGCCTTT");
#ifdef AGE_TIME
				timeval ali_s,ali_e;
				gettimeofday(&ali_s,NULL);
#endif
				Scorer scr(ops.match, ops.mismatch, ops.gapOpen, ops.gapExtend);
				if (ops.both) {
					Sequence *seq3 = seq2->clone();
					seq3->revcom();
					AGEaligner aligner1(*seq1,*seq2);
					AGEaligner aligner2(*seq1,*seq3);
					bool res1 = aligner1.align(scr,flag);
					bool res2 = aligner2.align(scr,flag);
					if (!res1 && !res2) { 
						cerr<<"No alignment made."<<endl;
					} else if (aligner1.score() >= aligner2.score()) {
						aligner1.printAlignment();
					} else {
						aligner2.printAlignment();
					}
					delete seq3;
				} else {
					AGEaligner aligner(*seq1,*seq2);
					if (aligner.align(scr,flag)) aligner.printAlignment();
					else cerr<<"No alignment made."<<endl;
				}
#ifdef AGE_TIME
				gettimeofday(&ali_e, NULL);
				cout << "\nAlignment time is " << ali_e.tv_sec - ali_s.tv_sec +
                        (ali_e.tv_usec - ali_s.tv_usec)/1e+6<<" s"<<endl<<endl;
#endif
			}
		}
		Sequence::deleteSequences(old_seq1);
		Sequence::deleteSequences(old_seq2);
	}
	cerr << "\n******************** AGE ALL DONE ********************\n";

	return 0;
}

