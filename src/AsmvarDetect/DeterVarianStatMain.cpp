/* Author : Shujia Huang
 * Date   : 2014-02-26
 *
 * The main part of variant detection program
 *
 */

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <set>
#include "utility.h"
#include "Variant.h"
#include "gzstream.h"

using namespace std;

void Usage         ( const char* prog );
void ReadFileList  ( const char* filelist,    vector< string   > & infile );
void LoadBedRegion ( const string infile, map< string, vector<Region> > & region, map< string, size_t> & index );
void LoadAlign     ( vector< string > infile, AxtVar & variant );
map< string, vector<Region> > VariantRegs ( const string infile);
map<string,map<size_t,set<string> > > Overlap (map<string,vector<Region> >& tarReg, map<string,vector<Region> >& varReg, map<string,size_t> indx, bool isC);

int main ( int argc, char* argv[] ) {

	char c;
    vector< string > infile, tarRef, qryRef;
    string filelist, bedInfile, varInfile;
    while ( (c = getopt( argc, argv, "i:l:v:b:t:q:h" )) != -1 ) {
        switch ( c ){
            case 'i' : infile.push_back(optarg); break; // lastz axt file
            case 't' : tarRef.push_back(optarg); break; // Target targeterence in axt result. [fa format]
            case 'q' : qryRef.push_back(optarg); break; // Query  targeterence in axt result. [fa format]
            case 'l' : filelist  = optarg;       break; // lastz axt file list
            case 'b' : bedInfile = optarg;       break;
            case 'v' : varInfile = optarg;       break;
            case 'h' : Usage( argv[0] );
            default  :
                cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl;
                Usage( argv[0] );
        }
    }
    cerr << "Parameter : " << join(argv, argc) << "\n" << endl;
    if ( !filelist.empty() ) ReadFileList( filelist.c_str(), infile );
    if ( argc == 1 || infile.empty() || bedInfile.empty() || varInfile.empty() || qryRef.empty() || tarRef.empty() ) Usage( argv[0] );

	map< string, vector<Region> > bedReg;
	map< string, size_t         > indx;
	LoadBedRegion(bedInfile, bedReg, indx);

    AxtVar variant;
    for ( size_t i (0); i < qryRef.size(); ++i ) { variant.qryfa.Load( qryRef[i] ); } // Load query sequence
    for ( size_t i (0); i < tarRef.size(); ++i ) { variant.tarfa.Load( tarRef[i] ); } // Load query sequence
	LoadAlign(infile, variant);
	//map< string, vector<Region> > varRegs = variant.VarTarRegs();
	map< string, vector<Region> > varRegs = VariantRegs( varInfile );

	map<string, map<size_t, set<string> > > overlap2mapTarNum = Overlap(variant.maptar, bedReg, indx, true);
	map<string, map<size_t, set<string> > > overlap2varTarNum = Overlap(varRegs, bedReg, indx, false);

	// Outputting
	for ( map< string,vector<Region> >::iterator it( bedReg.begin() ); it != bedReg.end(); ++it ) {
		string id = it->first;
		for (size_t i(0); i < bedReg[id].size(); ++i) {
			if ( overlap2mapTarNum.count(id) && overlap2mapTarNum[id].count(i) && overlap2varTarNum.count(id) && overlap2varTarNum[id].count(i) ) {
				for (set<string>::iterator is(overlap2varTarNum[id][i].begin()); is != overlap2varTarNum[id][i].end(); ++is) {
					if ( overlap2mapTarNum[id][i].count( *is ) && overlap2varTarNum[id][i].size() == 1 ) overlap2mapTarNum[id][i].erase(*is);
				}
			}
			int n1 = (overlap2mapTarNum.count(id) && overlap2mapTarNum[id].count(i) ) ? overlap2mapTarNum[id][i].size(): 0;
			int n2 = (overlap2varTarNum.count(id) && overlap2varTarNum[id].count(i) ) ? overlap2varTarNum[id][i].size(): 0;
			cout << id << "\t" << bedReg[id][i].start << "\t" << bedReg[id][i].end << "\t" << n1 << "\t" << n2 << "\t"
				 << bedReg[id][i].info << "\n";
		}
	}

	cerr << "*********************** The Program is Done ***********************" << endl;
	return 0;
}

map< string, vector<Region> > VariantRegs ( const string infile ) {

	map< string, vector<Region> > region;
	igzstream I ( infile.c_str() );
    Region reg;
	string tmp;
    while ( 1 ) {
	//20	16385	16670	286	SimulGap	C38846331	36	54
        I >> reg.id >> reg.start >> reg.end >> tmp >> tmp >> reg.info;;
        if ( I.eof() ){
            break;
        } else if ( I.fail() ) {
            cerr << "File read ERROR [I/O]!!\n" << reg.id << endl;
            exit(1);
        }
        getline (I, tmp, '\n');
        region[reg.id].push_back( reg );
    }
    I.close();

    for ( map< string,vector<Region> >::iterator it( region.begin() ); it != region.end(); ++it ) {
        sort ( region[it->first].begin(), region[it->first].end(), SortRegion );
    }
	
	return region;
}

map<string,map<size_t,set<string> > > Overlap ( map<string,vector<Region> >& tarReg, map<string,vector<Region> >& varReg, map<string, size_t> indx, bool isC ) {

	map<string, map<size_t, set<string> > > overlapNum;
	for ( map< string,vector<Region> >::iterator it( tarReg.begin() ); it != tarReg.end(); ++it ) {

		string id = it->first;
		if ( !indx.count(id) ) continue;
		for (size_t i(0); i < tarReg[id].size(); ++i) {
			
			bool flag(true);
			for (size_t j(indx[id]); j < varReg[id].size(); ++j) {
				if ( j > 0 ) assert( varReg[id][j-1].start <= varReg[id][j].start ); // Check sorted
				if ( tarReg[id][i].start > varReg[id][j].end  ) continue;
				if ( tarReg[id][i].end   < varReg[id][j].start) break;

				if ( flag ) {
					flag    = false;
					indx[id]= j;
				}

//if (varReg[id][j].start == 2377740) cerr << "# isC : " << isC << "\t" << id << " " << j << "\t" << tarReg[id][i].info << endl;
				if ( !isC ) {
					overlapNum[id][j].insert(tarReg[id][i].info);
				} else if ( tarReg[id][i].start <= varReg[id][j].start && tarReg[id][i].end >= varReg[id][j].end ) {
					overlapNum[id][j].insert(tarReg[id][i].info);
				}
			}
		}
	}

	return overlapNum;
}

void LoadBedRegion (const string infile, map< string, vector<Region> > & region, map< string, size_t> & index) {

	igzstream I ( infile.c_str() );
	Region reg;
	while ( 1 ) {
		I >> reg.id >> reg.start >> reg.end;
		if ( I.eof() ){
			break;
		} else if ( I.fail() ) {
			cerr << "File read ERROR [I/O]!!\n" << reg.id << endl;
			exit(1);
		}
		getline (I, reg.info, '\n');
		region[reg.id].push_back( reg );
	}
	I.close();

	for ( map< string,vector<Region> >::iterator it( region.begin() ); it != region.end(); ++it ) {
		sort ( region[it->first].begin(), region[it->first].end(), SortRegion );
		index[it->first] = 0;
	}
}

void LoadAlign (vector< string > infile, AxtVar & variant ) {

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

            variant.CheckAxt  ();
            variant.GetMapReg (); // The coordinate coversion events will be happen in this memerber function!
        }
        I.close();
    }

	for ( map< string,vector<Region> >::iterator it( variant.maptar.begin() ); it != variant.maptar.end(); ++it ) {
        sort ( variant.maptar[it->first].begin(), variant.maptar[it->first].end(), SortRegion );
    }
	return;
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

    cerr << "Version : 0.0.0 ( 2014-02-24 )                                                          \n"
        << "Author  : Shujia Huang                                                                   \n"
        << "Created : 2013-02-24                                                                     \n"
        << "\nCaution: The the input axt file should be the regular axt format. It means do not covert\n"
        << "           the coordinates of querys which align to '-' strand youself. It's not right and excrescent!\n"
        << "\nUsage : " << prog << " [Options] ""-i [axt infile] -o output                           \n"
        << "\n   Options  :                                                                          \n"
        << "       -i  [str]   axt file. require!                                                    \n"
        << "       -l  [str]   axt file list. [none]                                                 \n"
        << "       -b  [str]   Variant bed file. require!                                            \n"
        << "       -v  [str]   Variant bed file. require!                                            \n"
        << "       -t  [str]   target sequence, fa format. require!                                  \n"
        << "       -q  [str]   Query  sequence, fa format. require!                                  \n"
        << "       -h          Output this help information.                                         \n"
        << endl;

    exit(1);
}



