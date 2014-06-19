/* Author : Shujia Huang
 * Date   : 2013-12-30 10:12:53	Add gzstream to read file
 *          2013-10-18 15:25:18
 *
 * The main part of variant detection program
 *
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include "utility.h"
#include "Variant.h"
#include "gzstream.h"

using namespace std;

void Usage        ( const char* prog );
void ReadFileList ( const char* filelist,    vector< string   > & infile );

int main ( int argc, char* argv[] ) {

	char c;
	vector< string > infile, tarRef, qryRef;
	string filelist, outFilePrefix;
	while ( (c = getopt( argc, argv, "i:l:o:t:q:h" )) != -1 ) {
		switch ( c ){
            case 'i' : infile.push_back(optarg); break; // lastz axt file
            case 't' : tarRef.push_back(optarg); break; // Target targeterence in axt result. [fa format]
            case 'q' : qryRef.push_back(optarg); break; // Query targeterence in axt result. [fa format]
			case 'o' : outFilePrefix = optarg;   break;
            case 'l' : filelist      = optarg;   break; // lastz axt file list
            case 'h' : Usage( argv[0] );
            default  :
                cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl;
                Usage( argv[0] );
        }
	}
	cerr << "Parameter : " << join(argv, argc) << "\n" << endl;
	if ( !filelist.empty() ) ReadFileList( filelist.c_str(), infile );
	if ( argc == 1 || infile.empty() || outFilePrefix.empty() || qryRef.empty() || tarRef.empty() ) Usage( argv[0] );

	Variant variant;
	for ( size_t i (0); i < qryRef.size(); ++i ) { variant.qryfa.Load( qryRef[i] ); } // Load query sequence
	for ( size_t i (0); i < tarRef.size(); ++i ) { variant.tarfa.Load( tarRef[i] ); } // Load query sequence

	for ( size_t i (0); i < infile.size(); ++i ) {

		igzstream I ( infile[i].c_str() );
        if ( !I ){
            cerr << "Cannot open file : " << i + 1 << "\t" << infile[i] << endl;
            exit(1);
        }
        cerr << "***#    Reading file : " << i + 1 << "\t" << infile[i] << "    #*** " << endl;
		string tmp;
		while ( 1 ) {
		/*  1 chrM 16308 16389 Contig102837 1 81 + 5938
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
			if ( tmp[0] == '#' ){ getline ( I, tmp, '\n' ); continue; }
			variant.id.assign(tmp); // 1
			I >> variant.target.id >> variant.target.start >> variant.target.end; // chrM 16308 16389
			I >> variant.query.id  >> variant.query.start  >> variant.query.end ; // Contig102837 1 81
			I >> variant.strand    >> variant.score;                              // + 5938
			getline ( I, tmp, '\n' );
			getline ( I, variant.tarSeq, '\n' );  
			getline ( I, variant.qrySeq, '\n' );
			getline ( I, tmp, '\n' ); // Space line

			variant.CheckAxt();
			variant.CallSNP ();
			variant.CallInsertion(); // Don't covert the coordinate which map to the '-' strand here. I'll covert it when calling indel.
			variant.CallDeletion (); // Don't covert the coordinate which map to the '-' strand here. I'll covert it when calling indel.
			variant.GetMapReg    (); // The coordinate coversion events will be happen in this memerber function!
		}
		I.close();
	}
	variant.CallSV(); // It's SV not Indel, which cannot be called by a single alignment. The most important part of these program!!
	variant.CallClipReg();
	variant.CallNomadic();
	variant.Filter();        // Filter the indels' regions which in nosolution regions

	variant.Output   ( outFilePrefix + ".svd"     );
	variant.Summary  ( outFilePrefix + ".summary" );
	variant.OutputSNP( outFilePrefix + ".snp"     );
	variant.OutputGap( outFilePrefix + ".gap.bed" );

	cerr << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> All Done <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	return 0;
}

void ReadFileList ( const char* filelist, vector< string > & infile ) {

    igzstream I ( filelist );
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

void Usage ( const char* prog ) {

    cerr << "Version : 0.0.0 ( 2013-10-18 )                                                          \n"
        << "Author  : Shujia Huang                                                                   \n"
        << "Created : 2013-10-18                                                                     \n"
        << "Last Modify : 2013-11-19 17:31:49 Fix bugs                                               \n"
		<< "\nCaution: The the input axt file should be the regular axt format. It means do not covert\n"
		<< "           the coordinates of querys which align to '-' strand youself. It's not right and excrescent!\n"
        << "\nUsage : " << prog << " [Options] ""-i [axt infile] -o output                           \n"
        << "\n   Options  :                                                                          \n"
        << "       -i  [str]   axt file. require!                                                    \n"
        << "       -l  [str]   axt file list. [none]                                                 \n"
        << "       -o  [str]   Output file prefix. require!                                          \n"
        << "       -t  [str]   target sequence, fa format. require!                                  \n"
        << "       -q  [str]   Query  sequence, fa format. require!                                  \n"
        << "       -h          Output this help information.                                         \n"
        << endl;

    exit(1);
}

