/*
 *
 *  Author : Shujia Huang
 *  Date   : 2013/05/24
 *	
 *  Last modify : 
 *  	2013-10-2 	Update
 *  	2013-09-21	Add a GetDeletetion() to detect deletion.
		2013-07-25	Delete command parameter '-L', we can just use the query fa file to the length of query sequnce.
					Delete function. 'ReadQrylength()', cause we dontt need it any more.
					Add "tLengthList" to record the length of target.
					Add function 'GetFaLength()' to record the fa size.
		2013-6-26	Add output the totally unmapped query as 'nomadic'
		2013-05-31	Update the parameter '-t' and '-q', make them able to get multiple target fa files.
 *
 *  This program is used to call Novel sequence region and gap region against the targeterence.
 *	This program can just calculate the the depeth lower than 127. if greater than 255, it'll treat to 127.
 *
 *  [Input] We need three input files. 
        1. lastz axt file
        2. length of query, format: HUMuvjD320poolingDFAAPEI-1_Index1_Index5_C3243  1389
            e.g. /ifs1/HP/PROJECT/YH_FOSMID/huangshujia/data/YHFosmid.fa.list
        3. Reference's N regions, format: chr1    177418  227417 
            e.g. /ifs1/ST_EPI/USER/huangshujia/database/human_targeterence/hg19.N_region.txt
    [Output] Output two files, one is for query, the other is for target sequence.  The format of these two files are : fa like.
 *
 */
#include <iostream>
#include <fstream>
#include <assert.h>
#include <getopt.h>
#include <vector>
#include <map>
#include <set>
#include <cctype>
#include <stdlib.h>
#include <utility>
//#include "/ifs1/ST_EPI/USER/huangshujia/bin/cpp_bin/lib/mystl/utility.h"
#include "utility.h"
#include "AxtVar.h"

using namespace std;

typedef struct {
    string id;
    int start;
    int end;
} Region;

typedef struct {

    Region target;
    Region query;

    char strand;
    long score;

    string tarBases;
    string qryBases;

} AXT;

typedef struct {

    Region target;
    Region query;
    string type;
} R2R;

void Usage        ( const char* prog );
void ReadFileList ( const char* filelist,    vector< string   > & infile );
void GetFaLength  ( map<string, string> &fa, map< string, int > & lengthList );
void ReadFaSeq    ( int mode, vector<string> & file, map<string, string> & fa, bool isMaskeN = false );
void LoadData     ( int mode, vector<string> & infile, map<string, int> & qLengthList, map<string, string> & tarFa, map<string, string> & qryFa, bool isCov );
void MaskerTwoEnd ( AXT & axt, string & tarFa, string & qryFa );
void Coverage     ( AXT & axt, string & tarFa, string & qryFa ); // Calculate the coverage of axt mapping result.
void OutAxtErr    ( AXT & info );
void OutErrRegion ( vector< R2R > & region );
void Output       ( map< string, string > & fa, const char* outfile );
bool Check        ( vector<R2R> & region );
void GetRegionIntrAXT            ( AXT & info, unsigned int tarLen, string & qryFa );
void GetDeletion                 ( AXT & info, unsigned int qryLen, string & tarFa );
void GetNovleRegionAndClipRegion ( map< string, vector<R2R> > & mapInfo, map<string, string> & tarFa, map<string, int> & qLengthList );

int main( int argc, char* argv[] ) {

	char c;
	bool isCov( true );
	int  mode(1);
	vector< string > infile, tarRef, qryRef;
	string filelist, outfilePrefix /*, queryLengthListFile*/;
	while ( (c = getopt( argc, argv, "i:l:o:t:q:m:fh" )) != -1 ){
		switch ( c ){
			case 'i' : infile.push_back(optarg); break; // lastz axt file 
			case 't' : tarRef.push_back(optarg); break; // Target targeterence in axt result.[fa format]
			case 'q' : qryRef.push_back(optarg); break; // Query targeterence in axt result. [fa format]
			case 'l' : filelist      = optarg;   break; // lastz axt file list
			case 'o' : outfilePrefix = optarg;   break;
			case 'f' : isCov         = false;    break;
			case 'm' : mode          = atoi(optarg); break;
			//case 'L' : queryLengthListFile = optarg; break;
			case 'h' : Usage( argv[0] );
			default  : 
				cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl; 
				Usage( argv[0] );
		}
	}
	cerr << "Parameter : " << join(" ", argv, argc ) << endl;
	if ( !filelist.empty() ) ReadFileList( filelist.c_str(), infile );
	if ( argc == 1 || infile.empty() || (outfilePrefix.empty() && mode==1) || tarRef.empty() || qryRef.empty() ) Usage( argv[0] );


	map< string, string > tarFa, qryFa; // Target fa, Query fa.
	ReadFaSeq ( mode, tarRef, tarFa );
	ReadFaSeq ( mode, qryRef, qryFa, true );

	map< string, int > qLengthList, tLengthList; // Query length list, Target length list[In fact I don't use target lengthi].
	GetFaLength ( qryFa, qLengthList ); // Query fa length

	cerr << "\n**************************** Processing Data *************************" << endl;
	LoadData ( mode, infile, qLengthList, tarFa, qryFa, isCov ); // AXT format.

	cerr << "\n**************************** Outputing  Data for model 1 *************************" << endl;
	if ( mode == 1 ) { // Just used in mode 1
		string outfile1( outfilePrefix + ".query.depth" ), outfile2( outfilePrefix + ".target.depth" );
		Output( qryFa, outfile1.c_str() );
		Output( tarFa, outfile2.c_str() );
		cerr << "The final result : \n\t" << outfile1 << "\n\t" << outfile2 << endl;
	}
	cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> All  Done <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	return 0;
}

void LoadData ( int mode, vector< string > & infile, map< string, int > & qLengthList, map< string, string > & tarFa, map< string, string > & qryFa, bool isCov ) {

	assert ( mode == 1 || mode == 2 );

	R2R mpI;
	set< string              > mapQuery;
	map< string, vector<R2R> > mapInfo;

	AXT info;
	string tmp;
	for ( size_t i (0); i < infile.size(); ++i ) {

		ifstream I ( infile[i].c_str() );
		if ( !I ){
			cerr << "Cannot open file : " << i + 1 << "\t" << infile[i] << endl;
			exit(1);
		}
		cerr << "Reading file : " << i + 1 << "\t" << infile[i] << endl;

		while ( 1 ) {
			/*
			   1 chrM 16308 16389 Contig102837.1_6242-11697_5456 1 81 + 5938
			   CATAGTACATAAAGTCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCCTCAGATA
			   CATAGCACATATAGTCATTCATCGTACATAGCACATTATAGTCAAATCATTTCTCGTCCCCACGGAT-ATCCCCCTCAGATA

			 */
			I >> tmp;
			if ( I.eof() ){
				break;
			} else if ( I.fail() ) {
				cerr << "File read ERROR [I/O]!!\n" << tmp << endl;
				exit(1);
			}

			if ( tmp[0] == '#' ){

				getline ( I, tmp, '\n' );
				continue;
			}
			I >> info.target.id; I >> info.target.start; I >> info.target.end;
			I >> info.query.id;  I >> info.query.start;
			I >> info.query.end; I >> info.strand >> info.score;

			getline ( I, tmp, '\n' );
			getline ( I, info.tarBases, '\n' );
			getline ( I, info.qryBases, '\n' );
			getline ( I, tmp, '\n' ); // Space line

			mapQuery.insert( info.query.id );
			if ( !qLengthList.count( info.query.id ) || !qryFa.count( info.query.id ) ) {

				cerr << "Miss some query id or query id can't match!!!" << endl;
				cerr << "The unmatch query(main) : " << info.query.id   << endl;
				OutAxtErr( info ); exit (1);
			}
			if ( !tarFa.count( info.target.id ) ) {

				cerr << "Miss some target ids or target id can't match!!!" << endl;
				cerr << "The unmatch target: " << info.target.id           << endl;
				OutAxtErr( info ); exit (1);
			}

			if ( info.strand == '-' && isCov ) { //isCov is just used here

				int itmp = info.query.end;
				info.query.end   = qLengthList[info.query.id] - info.query.start + 1;
				info.query.start = qLengthList[info.query.id] - itmp + 1;
			}
			mpI.target.id = info.target.id; mpI.target.start = info.target.start; mpI.target.end = info.target.end;
			mpI.query.id  = info.query.id ; mpI.query.start  = info.query.start ; mpI.query.end  = info.query.end;
			mapInfo[mpI.target.id + "#" + mpI.query.id].push_back ( mpI );

			if ( mode == 1 ) {
				// Continue here ... "learn soap.coverage" output fa like. 
				// I have to calculate the fasta query sequence, so that I can get the novel sequence easily.
				if ( (qLengthList[info.query.id] == info.query.end) || (info.query.start == 1) ) { 
					MaskerTwoEnd(info,tarFa[info.target.id],qryFa[info.query.id]);
				}
				Coverage( info, tarFa[info.target.id], qryFa[info.query.id] );
			} else if ( mode == 2 ) {

				GetRegionIntrAXT ( info, tarFa[info.target.id].length(), qryFa[info.query.id] );
				GetDeletion ( info, qryFa[info.query.id].length(), tarFa[info.target.id] );
			} else {
				cerr << "Mode ERROR: " << mode << endl;
				exit(1);
			}
		}
		I.close();
	}

	cerr << endl;
	if ( mode == 2 ) GetNovleRegionAndClipRegion ( mapInfo, tarFa, qLengthList );

	// Output nomadic 
	for ( map<string, int>::iterator it ( qLengthList.begin() ); it != qLengthList.end(); ++it ) {

		if ( mapQuery.count( it->first ) ) continue;
		cout << "-\t-\t-\t-\t" << it->first         << "\t1\t"        << it->second << "\t" 
			 << it->second    << "\t" << it->second << "\tnomadic\tN" << endl; 
	}
	cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>    All files Done   <<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	return;
}

void GetRegionIntrAXT ( AXT & info, unsigned int tarLen, string & qryFa ) {

	vector< R2R > region;

    if ( info.tarBases.size() != info.qryBases.size() ){
		cerr << "Size Unmatch!!\n" << info.target.id << "\t" << info.target.start << "\t"  << info.query.id << "\t"
			 << info.query.start   << "\nRef_base: " << info.tarBases << "\nQuery_bases: " << info.qryBases
			 << "\nProgram break!" << endl;
        exit(1);
    }

    int tarPosIncreat (0);
    int qryPosIncreat (0);
    int tarPos ( info.target.start ), gapRefStart, gapRefEnd;
    int qryPos ( info.query.start  ), gapQryStart, gapQryEnd;

    bool findGapS( false ), findGapE( false );
    R2R r2r;
    for ( size_t i(0); i < info.tarBases.size(); ++i ) {

        tarPos += tarPosIncreat;
        qryPos += qryPosIncreat;
        tarPosIncreat = ( info.tarBases[i] == '-' ) ? 0 : 1;
        qryPosIncreat = ( info.qryBases[i] == '-' ) ? 0 : 1;

        if ( findGapS && ( (toupper(info.tarBases[i]) != 'N' && info.tarBases[i] != '-') || (info.qryBases[i] == '-') ) ) {
            findGapE  = true;
            gapRefEnd = tarPos; 
			gapQryEnd = qryPos;
        }
        if ( findGapS && findGapE ) {
            findGapS  = false;  findGapE  = false;
            r2r.target.id    = info.target.id; r2r.query.id   = info.query.id;
            r2r.target.start = gapRefStart;    r2r.target.end = gapRefEnd;
            r2r.query.start  = gapQryStart;    r2r.query.end  = gapQryEnd;
			r2r.type         = ".";
            if ( toupper(info.tarBases[i-1]) == 'N' ) r2r.type = "N";
            if ( info.tarBases[i-1] == '-' )          r2r.type = "-";
            region.push_back( r2r );
        }
//cerr << "# " << info.target.id << "\t" << tarPos << "\t" << info.tarBases[i] << " - " << info.qryBases[i] << "\t[" << qryFa[qryPos-1] << "]\t" << qryPos << "\t" << info.query.id << endl;

		if ( info.qryBases[i] == '-' ) continue;

		info.qryBases[i] = qryFa[qryPos-1];
        if ( toupper (info.tarBases[i]) == 'N' || info.tarBases[i] == '-' ) {
            if ( !findGapS && toupper(info.qryBases[i]) != 'N' ) {
                findGapS    = true;
                //gapRefStart = ( info.tarBases[i] == '-' ) ? tarPos - 1 : tarPos;
                gapRefStart = tarPos - 1;
                gapQryStart = qryPos - 1;
//cerr << "[Gap] Start: " << gapRefStart << "\tgapQryStart: " << gapQryStart << endl;
            } else if ( findGapS && toupper(info.qryBases[i]) == 'N' ) {
                findGapE  = true;
                gapRefEnd = ( info.tarBases[i] != '-' ) ? tarPos + 1 : tarPos;
                gapQryEnd = qryPos + 1;
//cerr << "[Gap] End: " << gapRefEnd << "\tgapQryEnd: " << gapQryEnd << endl;
            }
        }
    }

    if ( findGapS && (toupper(info.tarBases.rbegin()[0]) == 'N' || info.tarBases.rbegin()[0] == '-') ) {
        r2r.target.id    = info.target.id;  r2r.query.id   = info.query.id;
        r2r.target.start = gapRefStart;     r2r.target.end = tarPos;
        r2r.query.start  = gapQryStart;     r2r.query.end  = qryPos;
		r2r.type         = ".";
		if ( toupper(info.tarBases.rbegin()[0]) == 'N' ) r2r.type = "N";
		if ( info.tarBases.rbegin()[0] == '-' )          r2r.type = "-";
        region.push_back(r2r);
    }

	for ( size_t i(0); i < region.size(); ++i ) {

		cout << region[i].target.id << "\t" << region[i].target.start << "\t" << region[i].target.end << "\t" << tarLen << "\t"
			 << region[i].query.id  << "\t" << region[i].query.start  << "\t" << region[i].query.end  << "\t"
			 << region[i].query.end - region[i].query.start - 1       << "\t" << qryFa.length()       << "\t" << region[i].type << "\t" 
			 << info.strand << endl;
	}

    return;
}

void GetDeletion ( AXT & info, unsigned int qryLen, string & tarFa ) {

    vector< R2R > region;

    if ( info.tarBases.size() != info.qryBases.size() ){
        cerr << "Size Unmatch!!\n" << info.target.id << "\t" << info.target.start << "\t"  << info.query.id << "\t"
             << info.query.start   << "\nRef_base: " << info.tarBases << "\nQuery_bases: " << info.qryBases
             << "\nProgram break!" << endl;
        exit(1);
    }

    int tarPosIncreat (0);
    int qryPosIncreat (0);
    int tarPos ( info.target.start ), gapRefStart, gapRefEnd;
    int qryPos ( info.query.start  ), gapQryStart, gapQryEnd;

    string tarStr( info.tarBases ), qryStr( info.qryBases );

    bool findGapS( false ), findGapE( false );
    R2R r2r;
	for ( size_t i(0); i < tarStr.size(); ++i ) {

        tarPos += tarPosIncreat;
        qryPos += qryPosIncreat;
        tarPosIncreat = ( tarStr[i] == '-' ) ? 0 : 1;
        qryPosIncreat = ( qryStr[i] == '-' ) ? 0 : 1;

        if ( findGapS && ( (qryStr[i] != '-') || (tarStr[i] == '-') ) ) {
            findGapE  = true;
            gapRefEnd = tarPos;
            gapQryEnd = qryPos;
        }
        if ( findGapS && findGapE ) {
            findGapS  = false;  findGapE  = false;
            r2r.target.id    = info.target.id; r2r.query.id   = info.query.id;
            r2r.target.start = gapRefStart;    r2r.target.end = gapRefEnd;
            r2r.query.start  = gapQryStart;    r2r.query.end  = gapQryEnd;

            r2r.type = "D";
            region.push_back( r2r );
        }
        if ( tarStr[i] == '-' ) continue;

        tarStr[i] = tarFa[tarPos-1];
        if ( qryStr[i] == '-' ) {
            if ( !findGapS && toupper(tarStr[i]) != 'N' ) {
                findGapS    = true;

                gapRefStart = tarPos - 1;
                gapQryStart = qryPos - 1;

            } else if ( findGapS && toupper(tarStr[i]) == 'N' ) {
                findGapE  = true;
                gapQryEnd = ( qryStr[i] != '-' ) ? qryPos + 1 : qryPos;
                gapRefEnd = tarPos + 1;

            }
        }
    }

	if ( findGapS && (toupper(tarStr.rbegin()[0]) == 'N' || tarStr.rbegin()[0] == '-') ) {
        r2r.target.id    = info.target.id;  r2r.query.id   = info.query.id;
        r2r.target.start = gapRefStart;     r2r.target.end = tarPos;
        r2r.query.start  = gapQryStart;     r2r.query.end  = qryPos;

        r2r.type         = "D";
        region.push_back(r2r);
    }

    for ( size_t i(0); i < region.size(); ++i ) {

        cout << region[i].target.id << "\t" << region[i].target.start << "\t" << region[i].target.end << "\t" << tarFa.length() << "\t"
             << region[i].query.id  << "\t" << region[i].query.start  << "\t" << region[i].query.end  << "\t"
             << region[i].query.end - region[i].query.start - 1       << "\t" << qryLen               << "\t" << region[i].type << "\t"
             << info.strand << endl;
    }

    return;
}

void GetNovleRegionAndClipRegion ( map< string, vector<R2R> > & mapInfo, map<string, string> & tarFa, map<string, int> & qLengthList ) {


	map < string, vector<R2R> > qryMapInfo; // Just for getting clip region
	for ( map< string, vector<R2R> >::iterator it( mapInfo.begin() ); it != mapInfo.end(); ++it ) {

		sort( it->second.begin(), it->second.end(), mySortByTar );

		vector<R2R> unmapRegion;
		int tLeftStart, tLeftEnd, tRightStart, tRightEnd;
		int qLeftStart, qLeftEnd, qRightStart, qRightEnd;

		if ( !Check(it->second) ) { 
			OutErrRegion ( it->second );
			cerr << endl;

			for ( size_t i(0); i < it->second.size(); ++i ) { qryMapInfo[it->second[i].query.id].push_back( it->second[i] ); }
			continue;
		}
		// Get novel haplotype region
		for ( size_t i(0); i < it->second.size(); ++i ) {
			if ( !qLengthList.count(it->second[i].query.id) ) {
				cerr << "[ERROR]Miss some query id or query id can't match!!! Check your query list!" << endl;
				cerr << "The unmatch query : " << it->second[i].query.id.length() << "\t" << it->second[i].query.id << endl;
				exit(1);
			}
			qryMapInfo[it->second[i].query.id].push_back( it->second[i] ); // Use for getting query mapped region to get clip region. 

			tRightStart = it->second[i].target.start; tRightEnd = it->second[i].target.end;
			qRightStart = it->second[i].query.start ; qRightEnd = it->second[i].query.end ;

			if ( i > 0 ) {

				if ( qLeftStart > qRightStart ) {
				
					cout<< it->second[i].target.id << "\t" << tLeftEnd  << "\t" << tRightStart    << "\t" 
					    << tarFa[it->second[i].target.id].length()                                << "\t"
						<< it->second[i].query.id  << "\t" << qRightEnd << "\t" << qLeftStart     << "\t"
						<< qLeftStart - qRightEnd-1<< "\t" << qLengthList[it->second[i].query.id] << "\tnovel\tN" << endl; // includeing novel haplotype and insertions.
				} else {

					cout<< it->second[i].target.id << "\t" << tLeftEnd << "\t" << tRightStart     << "\t"
						<< tarFa[it->second[i].target.id].length()                                << "\t"
						<< it->second[i].query.id  << "\t" << qLeftEnd << "\t" << qRightStart     << "\t"
						<< qRightStart - qLeftEnd-1<< "\t" << qLengthList[it->second[i].query.id] << "\tnovel\tN" << endl; // includeing novel haplotype and insertions.
				}
			}
			tLeftStart = tRightStart; tLeftEnd = tRightEnd;
			qLeftStart = qRightStart; qLeftEnd = qRightEnd;
		}
	}

	// Getting clip region
	for ( map< string, vector<R2R> >::iterator it( qryMapInfo.begin() ); it != qryMapInfo.end(); ++it ) {

		sort( it->second.begin(), it->second.end(), mySortByQry );

		if ( it->second.front().query.start > 1 )
			cout << "-\t-\t-\t-\t" << it->second.front().query.id << "\t1\t"  << it->second.front().query.start - 1 << "\t"
                 << it->second.front().query.start - 1 << "\t"    << qLengthList[it->first]       << "\tclip\tN"    << endl;
		
		if ( it->second.back().query.end < qLengthList[it->first] )
			cout << "-\t-\t-\t-\t" << it->second.front().query.id << "\t" << it->second.back().query.end + 1
                 << "\t"           << qLengthList[it->first]      << "\t" << qLengthList[it->first] - it->second.back().query.end
				 << "\t"           << qLengthList[it->first]
                 << "\tclip\tN"    << endl;
	}
}
// For mode 2
bool Check ( vector<R2R> & region ) {

	if ( region.size() < 3 ) return true;

	bool flag( true ), cmp( false );

	int qPreStart = region[0].query.start;

	for ( size_t i(1); i < region.size() - 1; ++i ) {

		if ( (region[i].query.start <= qPreStart) && (region[i].query.start <= region[i+1].query.start) ) {
			flag = false;
		} else if ( (region[i].query.start >= qPreStart) && (region[i].query.start >= region[i+1].query.start) ) {
			flag = false;
		}
		qPreStart = region[i].query.start;
		if ( !flag ) break;
	}
	return flag;
}

// For mode 1 
void MaskerTwoEnd ( AXT & axt, string & tarFa, string & qryFa ) {

	if ( axt.tarBases.length() != axt.qryBases.length() ) {
		cerr << "target and query length is not matching. target: " << axt.tarBases.length() << "\tQuery: " << axt.qryBases.length() << endl;
        OutAxtErr( axt ); exit(1);
	}
	bool flag (false);
	if ( axt.tarBases[0] == '-' || axt.tarBases[0] == 'n' || axt.tarBases[0] == 'N' ) flag = true;
	if ( axt.tarBases[axt.tarBases.size()-1] == '-' || axt.tarBases[axt.tarBases.size()-1] == 'n' || axt.tarBases[axt.tarBases.size()-1] == 'N') flag = true;
	if ( flag ) {
		cerr << "\t*** Masking end event @: " << axt.target.id   << "\t" << axt.target.start << "\t" << axt.target.end << "\t"
			 << axt.query.id  << "\t"         << axt.query.start << "\t" << axt.query.end    << "\t" << axt.strand     << "\t" << axt.score
			 << endl;
	}

	int tarPos, qryPos, tarPosInc(0), qryPosInc(0);	

	// Masker header
	if ( axt.tarBases[0] == '-' || axt.tarBases[0] == 'n' || axt.tarBases[0] == 'N' ) {
		tarPos = axt.target.start; 
		qryPos = axt.query.start;

		for ( size_t i(0); i < axt.tarBases.length() && (axt.tarBases[i] == '-' || axt.tarBases[i] == 'n' || axt.tarBases[i] == 'N'); ++i ){
			tarPos += tarPosInc; tarPosInc = ( axt.tarBases[i] == '-' ) ? 0 : 1;
			qryPos += qryPosInc; qryPosInc = ( axt.qryBases[i] == '-' ) ? 0 : 1;

			if ( axt.tarBases[i] == '-' ) qryFa[qryPos-1]  = '-';
			if ( axt.tarBases[i] == 'n' || axt.tarBases[i] == 'N' ) tarFa[tarPos-1] = '-';
		}
	}

	// Masker tail
	if ( axt.tarBases[axt.tarBases.size()-1] == '-' || axt.tarBases[axt.tarBases.size()-1] == 'n' || axt.tarBases[axt.tarBases.size()-1] == 'N' ) {
		tarPos = axt.target.end;
		qryPos = axt.query.end;
		
		tarPosInc = 0; qryPosInc = 0;
		for ( size_t i(axt.tarBases.size() - 1); i >= 0 && (axt.tarBases[i] == '-' || axt.tarBases[i] == 'n' || axt.tarBases[i] == 'N'); --i ){
			tarPos += tarPosInc; tarPosInc = ( axt.tarBases[i] == '-' ) ? 0 : -1;
            qryPos += qryPosInc; qryPosInc = ( axt.qryBases[i] == '-' ) ? 0 : -1;

			if ( axt.tarBases[i] == '-' ) qryFa[qryPos-1]  = '-';
			if ( axt.tarBases[i] == 'n' || axt.tarBases[i] == 'N' ) tarFa[tarPos-1] = '-';
		}
	}
	return;
}

void Coverage ( AXT & axt, string & tarFa, string & qryFa ) {

	const int overflowTrd (127); // Overflow threshold.
	if ( axt.tarBases.length() != axt.qryBases.length() ) {
		cerr << "target and query length is not matching. target: " << axt.tarBases.length() << "\tQuery: " << axt.qryBases.length() << endl;
		OutAxtErr( axt ); exit(1);
	}

	int tarPosInc (0); // target position increase.
	int qryPosInc (0); // query  position increase.
	int tarPos( axt.target.start ), qryPos( axt.query.start );

	for ( size_t i(0); i < axt.tarBases.length(); ++i ) {

		tarPos += tarPosInc;
		qryPos += qryPosInc;
		tarPosInc = ( axt.tarBases[i] == '-' ) ? 0 : 1;
		qryPosInc = ( axt.qryBases[i] == '-' ) ? 0 : 1;

if ( qryPos > qryFa.length() ) qryPos = qryFa.length();

		if ( axt.tarBases[i] != '-' && axt.qryBases[i] != '-' ) { //Do not include refernce Novel(deletion) , gap regions and unmapped regions
			if ( tarPos > tarFa.length() || qryPos > qryFa.length() ) {
				cerr << "[ERROR]Why is the length of fa shorter than the position?! I think, your input fa may be wrong!!" << endl;
				cerr << "Target length: " << tarFa.length() << ". While position: " << tarPos << endl;
				cerr << "Query  length: " << qryFa.length() << ". While position: " << qryPos << endl;
				OutAxtErr( axt ); exit(1);
			}
			if ( (tarFa[tarPos-1] < overflowTrd) && (tarFa[tarPos-1] != '-') ) ++tarFa[tarPos-1];
			//if ( (qryFa[qryPos-1] < overflowTrd) && (qryFa[qryPos-1] != '-') ) ++qryFa[qryPos-1];
			if ( (qryFa[qryPos-1] < overflowTrd) && (qryFa[qryPos-1] != '-') ) { 
				if (tarFa[tarPos-1] != 'n' && tarFa[tarPos-1] != 'N') ++qryFa[qryPos-1];
			}
		}
	}

	return;
}

void ReadFileList ( const char* filelist, vector< string > & infile ) {

	ifstream I ( filelist );
	if ( !I ) { cerr << "Cannot open file : " << filelist << endl; exit(1); }
	string tmp;
	while ( 1 ) {

		I >> tmp;
		if ( I.eof() ) break;

		infile.push_back( tmp );
		getline( I, tmp, '\n');
	}
	I.close();

	return;
}

/*
void ReadQrylength ( const char* qLengthFile, map< string, int > &qLengthList ) {

	ifstream I( qLengthFile );
	if ( !I ) {
		cerr << "Cannot open file : " << qLengthFile << endl;
		exit(1);
	}
	qLengthList.clear();
	//HUMuvjD60poolingDDAAPEI-12_Index126_C27233_7180002912139   101    ......
	string ctg, lengthStr, tmp;
	while ( 1 ) {

		I >> ctg >> lengthStr;
		getline ( I, tmp, '\n' );
		if ( I.eof() ) {
			break;
		} else if ( I.fail() ) {
			cerr << "I/O ERROR" << endl;
		}
		qLengthList[ctg] = atoi ( lengthStr.c_str() );
	}
	I.close();
}
*/
void GetFaLength ( map<string, string> &fa, map< string, int > &lengthList ) {

	lengthList.clear();
	for ( map<string, string>::iterator it(fa.begin()); it != fa.end(); ++it ) {
		lengthList[it->first] = it->second.length();
	}
	return;
}

void ReadFaSeq ( int mode, vector< string > & file, map<string, string> & fa, bool isMaskeN ) {

	if ( !isMaskeN ) { cerr << "Do not maske the N region in : " << join ("\n", file) << endl; }

	for ( size_t i(0); i < file.size(); ++i ) {

		ifstream I( file[i].c_str() );
		if ( !I ){
			cerr << "Cannot open file : " << file[i].c_str() << endl;
			exit(1);
		}

		string tmp, faId;
		while ( 1 ) {

			I >> tmp;
			if ( I.eof() ) break;
			if ( tmp[0] != '>' ) {

				fa[faId] += tmp;
			} else {
				faId.assign ( tmp, 1, string::npos );
				fa[faId].clear();
			}
			getline ( I, tmp, '\n' );
		}
		I.close();
	}

	if ( mode == 1 ) {
		for ( map<string, string>::iterator it( fa.begin() ); it != fa.end(); ++it ) {

			if ( !isMaskeN ) {
				// assign the characters to 1, instead of 0, cause if 0,it'll be ambiguous in (find_first_not_of(0) in Output() function.).
				it->second.assign( it->second.length(), 1 );
			} else {
				for ( size_t i(0); i < it->second.size(); ++i ) {
					if ( it->second[i] != 'n' && it->second[i] != 'N' ) {
						// assign the characters to 1, instead of 0, cause if 0,it'll be ambiguous in (find_first_not_of(0) in Output() function.).
						it->second[i] = 1;
					} else {
						it->second[i] = '-'; //Maske N region.
					}
				}
			}
		}
	}

	return;
}

void Output ( map< string, string > & fa, const char* outfile ) {

	ofstream O( outfile );
	if ( !O ){ cerr << "Cannot write to file: " << outfile << endl; }
	for ( map<string, string>::iterator it( fa.begin() ); it != fa.end(); ++it ) {

//		if ( it->second.find_first_not_of( 1 ) == string::npos ) continue; // Do not output the fa which don't be covered to save the disk space!

		O << ">" << it->first << endl; // faId
		for ( size_t i(0); i < it->second.length(); ++i ) {
			if ( it->second[i] != '-' ) {
				O << int( it->second[i] ) - 1 << " ";
			} else {
				O << it->second[i] << " ";
			}
			if ( (i + 1) % 100 == 0 ) O << endl;
		}
		if ( it->second.length() % 100 > 0 ) O << endl;
	}
	O.close();
}

void OutAxtErr ( AXT & info ) {

    cerr << "#" << info.target.id << "\t" << info.target.start << "\t" << info.target.end << "\t"
         << info.query.id         << "\t" << info.query.start  << "\t" << info.query.end  << "\n"
         << info.tarBases         << "\n" << info.qryBases     << endl;
    return;
}

void OutErrRegion ( vector< R2R > & region ) {

	for ( size_t i(0); i < region.size(); ++i ) {


		cerr << "@ "<< region[i].target.id << "\t" << region[i].target.start << "\t" << region[i].target.end << "\t"
			 << region[i].query.id         << "\t" << region[i].query.start  << "\t" << region[i].query.end  << endl;
	}
}

void Usage ( const char* prog ) {

	cerr << "Version : 0.0.1 ( 2013-05-24 )                                                          \n"
		<< "Author  : Shujia Huang                                                                   \n"
		<< "Created : 2013/05/24                                                                     \n"
		<< "Last Modify :                                                                            \n"
		<< "Usage : " << prog << " [Options] ""-i [axt infile] > output                              \n"
		<< "   Options  :                                                                            \n"
		<< "       -i  [str]   axt file. require!                                                    \n"
		<< "       -l  [str]   axt file list. [none]                                                 \n"
		<< "       -o  [str]   Output file prefix. require!                                          \n"
		<< "       -t  [str]   Target sequence, fa format. require!                                  \n"
		<< "       -q  [str]   Query  sequence, fa format. require!                                  \n"
		<< "       -f  [bool]  Convert the query position when mapping reverse chain. [true]: covert! \n"
		<< "       -m  [int]   Mode. Just mode 1 and mode 2. [1]                                     \n"
		<< "       -h          Output this help information.                                         \n"
		<< endl;

	exit(1);
}








