/* Author : Shujia Huang
 * Date   : 2013-10-18 15:25:18
 *
 * The main part of variant detection program
 *
 */
#include <iostream>
#include <fstream>
#include <getopt.h>
#include "Variant.h"
#include "utility.h"

using namespace std;

void Usage ( const char* prog );
void ReadFileList ( const char* filelist,    vector< string   > & infile );

int main ( int argc, char* argv[] ) {

	char c;
    bool isCov( true );
    vector< string > infile, qryRef;
    string filelist;
    while ( (c = getopt( argc, argv, "i:l:o:q:fh" )) != -1 ) {
        switch ( c ){
            case 'i' : infile.push_back(optarg); break; // lastz axt file
            case 'q' : qryRef.push_back(optarg); break; // Query targeterence in axt result. [fa format]
            case 'l' : filelist      = optarg;   break; // lastz axt file list
            case 'f' : isCov         = false;    break;
            case 'h' : Usage( argv[0] );
            default  :
                cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl;
                Usage( argv[0] );
        }
    }
    cerr << "Parameter : " << join(" ", argv, argc ) << endl;
    if ( !filelist.empty() ) ReadFileList( filelist.c_str(), infile );
    if ( argc == 1 || infile.empty() || qryRef.empty() ) Usage( argv[0] );

	AxtVar variant;
    for ( size_t i (0); i < qryRef.size(); ++i ) { variant.qryfa.Load( qryRef[i] ); } // Load query sequence
    for ( size_t i (0); i < infile.size(); ++i ) {

		ifstream I ( infile[i].c_str() );
        if ( !I ){
            cerr << "Cannot open file : " << i + 1 << "\t" << infile[i] << endl;
            exit(1);
        }
        cerr << "Reading file : " << i + 1 << "\t" << infile[i] << endl;
        string tmp;
	
		while ( 1 ) {
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
			if ( isCov ) variant.ConvQryCoordinate();

			variant.GetMapReg();
            variant.CallInsertion();
            variant.CallDeletion ();
		}
	}

	variant.Filter();
	variant.Output();
	return 0;
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

void Usage ( const char* prog ) {

    cerr << "Version : 0.0.1 ( 2013-10-18 )                                                          \n"
        << "Author  : Shujia Huang                                                                   \n"
        << "Created : 2013-10-18                                                                     \n"
        << "Last Modify :                                                                            \n"
        << "Usage : " << prog << " [Options] ""-i [axt infile] > output                              \n"
        << "   Options  :                                                                            \n"
        << "       -i  [str]   axt file. require!                                                    \n"
        << "       -l  [str]   axt file list. [none]                                                 \n"
        << "       -q  [str]   Query  sequence, fa format. require!                                  \n"
        << "       -f  [bool]  Convert the query position when mapping reverse chain. [true]: covert!\n"
        << "       -h          Output this help information.                                         \n"
        << endl;

    exit(1);
}

