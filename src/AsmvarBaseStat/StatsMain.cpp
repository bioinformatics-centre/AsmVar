/*
 * Author : Shujia Huang
 * Date   : 2013-11-22 14:37:13
 *
 * Modfiy from '~/Bin/software_pip/Reapr/Reapr_1.0.16/src/task_stats.cpp' wrote by Matin Hunt
 * 	2013-11-26 23:02:07	I shift to calculate the read coverage instead of fragment coverage!!
 *
 */
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <cstring>
#include <string>
#include <list>
#include <getopt.h>
#include <algorithm>

#include "Stats.h"

struct CmdLineOptions {

	uint16_t minMapQuality;  // ignore reads mapped with quality less than this
	uint32_t minAlignScore; 
	uint32_t iniReadLength;  // Just use to init the size of coverage deque
    string bamInfile;
	string gapInfile;
	string refGCfile;
	string gc2covfile;
};
// deals with command line options: fills the options struct
void ParseOptions(int argc, char** argv, CmdLineOptions& ops);
void Usage ( const char* prog );
  
int main ( int argc, char** argv ) {

	CmdLineOptions ops;
	ParseOptions(argc, argv, ops);
	cerr << "\nCommand Important Parameters: -q " << ops.minMapQuality << " -a " << ops.minAlignScore << " -l " << ops.iniReadLength << "\n" << endl;

	Stats stats;
	map<string, list<pair<long,long> > > refgaps;

	LoadGaps            (ops.gapInfile,       refgaps);
	//GetInsersizeFromBam (ops.bamInfile , stats.insert);
	//Tabix tigc(ops.refGCfile);
	//LoadGC2Cov(ops.gc2covfile, stats.gc2cov);

	// Go through input bam file getting local stats
	cerr << "************ Go through input bam file getting local stats **********\n" << endl;
	BamReader    bamReader;
    SamHeader    header;
    RefVector    references;

	if (!bamReader.Open(ops.bamInfile)) {
        cerr << ERROR_PREFIX << "Error opening bam file '" << ops.bamInfile << "'" << endl;
        return 1; 
    }
	header     = bamReader.GetHeader();
    references = bamReader.GetReferenceData();

	bool     firstRecord (true);
	int32_t  currentRefID(-1);
	//int32_t  alignScore;
	uint32_t  alignScore;
	//int32_t  isToosmall, isToobig;
    string   currentRefIDstring, rg, sa;
	BamAlignment al;

	unsigned long num = 0;
	//cout << "#CHROM\tPOS\tTOTAL_COV\tPROPER_PAIR_COV\tLOW_ALIGN_SCORE_PROPER_PAIR\tBAD_INSERT_COV\tWRONG_ORIETATION_COV\tSINGLE_END_COV\t"
    //     << "CLIP_AND_SA_COV\tEXPECTED_COV(Correct by GC contain)\n";
	cout << "#CHROM\tPOS\tTOTAL_COV\tPROPER_PAIR_COV\tLOW_ALIGN_SCORE_PROPER_PAIR\tCROSS_READ_COV\tBAD_INSERT_COV\tWRONG_ORIETATION_COV\tSINGLE_END_COV\t"
         << "CLIP_AND_SA_COV\n";
	while ( bamReader.GetNextAlignment(al) ) {

		++num; if ( num % 1000000 == 0 ) cerr << ">>> Reading lines number : " << num << endl;
		if ( !al.IsMapped() || al.IsDuplicate() || al.MapQuality < ops.minMapQuality ) continue;
//al.GetTag("AS", alignScore);
//cerr << "alignScore: " << alignScore << "\t# " << al.GetTag("AS", alignScore) << "\t# " << al.HasTag("AS") << endl;

		if (!al.GetTag("AS", alignScore)) { cerr << "Read " << al.Name << " doesn't have an alignment score AS:... Cannot continue" << endl; exit(1); }
		if (!al.GetTag("RG", rg)) { cerr << "Read " << al.Name << " doesn't have an alignment score RG:... Cannot continue "        << endl; exit(1); }
		//if (!stats.insert.count(rg)) {
		//	cerr << "Read " << al.Name << " doesn't have a RG id( " << rg << " ) at bam header ... Cannot continue " << endl; exit(1);
		//}

		// Deal with the case when we find a new reference sequence in the bam
		if ( currentRefID != al.RefID ) {

			if ( firstRecord ) {
				firstRecord = false;
			} else {
				OutputStats( currentRefIDstring, references[currentRefID].RefLength + 1, stats, refgaps );
			}

			currentRefID       = al.RefID;
			currentRefIDstring = references[al.RefID].RefName;

			unsigned long startPos = 1;
			stats.properPair       = Coverage( ops.iniReadLength, startPos );
			stats.pLowAlignScore   = Coverage( ops.iniReadLength, startPos );
			stats.crossSignal      = Coverage( ops.iniReadLength, startPos ); // to be a cross signal reads then the reads must be a proper-pair first!
			stats.badInserSizePair = Coverage( ops.iniReadLength, startPos );
			stats.wrongOrientation = Coverage( ops.iniReadLength, startPos );
			stats.single           = Coverage( ops.iniReadLength, startPos );
			stats.clipAlign        = Coverage( ops.iniReadLength, startPos );
			//LoadGC(tigc, currentRefIDstring, stats.refGC);
		}
		// print all stats to left 100bp of current read mapped position
		// It's very important to keep away from the current mapped position 100bp, this distance is for reads clip situation!
		OutputStats ( currentRefIDstring, al.Position - 100, stats, refgaps ); // Acturally I dn't have to set -100 now! 2013-11-26 22:57:32 
		// Now the coordinate in Coverage should equal to al.Position

		//isToosmall = stats.insert[rg].first - 3 * stats.insert[rg].second; if ( isToosmall < 0 ) isToosmall = 1;
        //isToobig   = stats.insert[rg].first + 3 * stats.insert[rg].second; 

		// Caution: Bamtools give the 0-base coordinate system which interval is half closed '(]', I shift it to 1-base
		// coordinate system which interval is closed '[]' below.
		// al.Position is 0-base; 
		// al.GetEndPosition() is 1-base by default
		if ( al.CigarData[0].Type     == 'S' && (al.Position+1) > 1 ) stats.clipAlign.Add((al.Position+1) - 1);
		if ( al.CigarData.back().Type == 'S' && al.GetEndPosition()+1 <= references[al.RefID].RefLength ) stats.clipAlign.Add(al.GetEndPosition() + 1);

		short pairOrientation = GetPairOrientation( al );

		if ( !al.IsMateMapped() || pairOrientation == UNPAIRED || pairOrientation == DIFF_CHROM ) {

			stats.single.Add (al.Position + 1, al.GetEndPosition()); // al.Position is 0-base, shift to 1-base
		} else if ( pairOrientation == INNER ) { //Correct orientation but the insertsize could be good or bad

			if ( al.IsProperPair() ) {

				if ( alignScore >= ops.minAlignScore ) { // stats.properPair

					// Should consider cross signal first!
					map< char, vector<pair<int32_t,int32_t> > > indelSignal = IndelSignalRefPos(al.Position + 1, al.CigarData);
					if ( indelSignal.count('I') ) {
						for ( size_t i(0); i < indelSignal['I'].size(); ++i ) 
							stats.crossSignal.Add( indelSignal['I'][i].first, indelSignal['I'][i].second );
					}
					if ( indelSignal.count('D') ) {
						for ( size_t i(0); i < indelSignal['D'].size(); ++i ) 
                            stats.crossSignal.Add( indelSignal['D'][i].first, indelSignal['D'][i].second );
					}
					if ( !indelSignal.count('I') && !indelSignal.count('D') ) stats.properPair.Add ( al.Position + 1, al.GetEndPosition() );
				} else {                                 // stats.pLowAlignScore
					stats.pLowAlignScore.Add( al.Position + 1, al.GetEndPosition() ); // read coverage
				}
			} else {
				stats.badInserSizePair.Add ( al.Position + 1, al.GetEndPosition() );
			}
		} else if ( pairOrientation == SAME || pairOrientation == OUTTER  ) { // Wrong Orientation

			stats.wrongOrientation.Add ( (al.Position + 1), al.GetEndPosition() );
		} else {
			cerr << ERROR_PREFIX << "Didn't expect this to happen... " << al.Name << endl; exit(1);
		}
	}
	bamReader.Close();
	cerr << ">>> Total lines number in bamInfile : " << num << endl;
	// print the remaining stats from the last ref sequence in the bam
	OutputStats( currentRefIDstring, references[currentRefID].RefLength + 1, stats, refgaps );
	cerr << "***************************** All Done ******************************\n" << endl;

	return 0;
}

void ParseOptions(int argc, char** argv, CmdLineOptions& ops) {

	// Set Default
	ops.minMapQuality = 30;
	ops.minAlignScore = 90;
	ops.iniReadLength = 1000;
	char c;
	while ( (c = getopt(argc, argv, "b:a:l:q:r:g:c:h") ) != -1 ) {

		switch ( c ) {
			case 'b' : ops.bamInfile     = optarg;       break;
			case 'r' : ops.refGCfile     = optarg;       break;
			case 'g' : ops.gapInfile     = optarg;       break;
			case 'c' : ops.gc2covfile    = optarg;       break; // gc contain to coverage
			case 'q' : ops.minMapQuality = atoi(optarg); break;
			case 'a' : ops.minAlignScore = atoi(optarg); break;
			case 'l' : ops.iniReadLength = atoi(optarg); break;
			case 'h' : Usage( argv[0] );
			default  : 
				cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl; exit(1);
		}
	}
	//if ( ops.bamInfile.empty() || ops.refGCfile.empty() || ops.gc2covfile.empty() || ops.gapInfile.empty() ) Usage ( argv[0] );
	if ( ops.bamInfile.empty() || ops.gapInfile.empty() ) Usage ( argv[0] );
}

void Usage ( const char* prog ) {

	cerr << "\nVersion : 0.0.0 (2013-11-25 22:33:53)                                 \n"
		 << "Author  : Shujia Huang                                                \n\n"
		 << "\nUsage : " << prog << " [Option] [-b bamInfile] > Output             \n\n"
		 << "     Options :                                                        \n\n"
		 << "     -b  bam input file. Required!                                      \n"
		 //<< "     -r  File of reference's GC contain. Required!                    \n"
		 << "     -g  File of reference gap' regions. Required!                      \n"
		 //<< "     -c  File of GC contain vs coverage. Required!                    \n"
		 << "     -q  Min Mapped quality.  [30]                                      \n"
		 << "     -a  Min alignment score. [90]                                      \n"
		 << "     -l  Init fragment length.[1000] Not Suggest to change this parameter.\n"
		 << "     -h  Show this help.                                                \n" << endl;

	exit(1);
}






