/*
 * Author : Shujia Huang
 * Date   : 2013-11-22 14:37:13
 *
 * Modfiy from '~/Bin/software_pip/Reapr/Reapr_1.0.16/src/task_stats.cpp' wrote by Matin Hunt
 *
 */

#include <iostream>
#include <fstream>
#include <utility>      // std::pair
#include <map>
#include <math.h>       /* sqrt */

#include "api/BamMultiReader.h"
#include "api/BamReader.h"
#include "utility.h"

using namespace std;
using namespace BamTools;

// insert size calculate 
bool GetInsersizeFromBam ( string bamInfile, map< string, pair<int, int> >& is, int n=20000 );// RG => pair<int, int>. 'first' is mean, 'second' is SD
pair<int, int> MeanAndStddev ( vector<int> & v );

int main ( int argc, char** argv ) {

	map< string, pair<int, int> > is;
	GetInsersizeFromBam ( argv[1], is, 200000 );
	return 0;
}

bool GetInsersizeFromBam ( string bamInfile, map< string, pair<int, int> >& is, int n ) {

	bool flag ( false );

	BamReader bamReader;
	if (!bamReader.Open(bamInfile)) {
        cerr << "Error opening bam file '" << bamInfile << "'" << endl;
        exit( 1 );
    }
	const SamHeader header     = bamReader.GetHeader();
	const RefVector references = bamReader.GetReferenceData();
	for (size_t i(0); i < references.size(); ++i ) { cerr << "# " << i << "\t" << references[i].RefName << "\t" << references[i].RefLength << endl; }

	if ( header.ReadGroups.IsEmpty() ) { 
		cerr << "Error in your bamfile " << bamInfile << " does not contain Read Group '@RG' on header." << endl;
		exit( 1 ); 
	}

	map<string, int> LN; // Recorde the number to calculate the insertsize 
	map<string, vector<int> > ins;
	
	unsigned long int toolong(0);
	for ( SamReadGroupConstIterator rgit = header.ReadGroups.Begin(); rgit != header.ReadGroups.End(); ++rgit ) {
		LN[rgit->ID] = n;
		toolong     += n;
	}
	toolong *= (10 * header.ReadGroups.Size()); // magnify 10 times of the @RG number. Just use to monitor error

	cerr << "***************** Calculating the insert size of RG *****************\n" << endl;

// 19  8002 
	BamAlignment al;
//	cerr << bamReader.HasIndex() << "\nJump to 19, 8002 : " << bamReader.Jump( 19, 8013 ) << endl;
//	if ( bamReader.Jump( 19, 8013 ) ) cerr << "Jump to 19, 8002 " << endl;

	while ( bamReader.GetNextAlignment(al) ) {

cerr << "Position: " << al.Name << "\t" << al.Position << "\t" << al.GetEndPosition() << "\tMatePosition: " << al.MatePosition << "\tCigae: " << al.CigarData[0].Type << " " << al.CigarData[0].Length << "\t" << al.RefID << "\t" << references[al.RefID].RefName << "\t" << references[al.RefID].RefLength << endl; continue;

		bool isEnd ( true );
		for ( map<string, int>::iterator it(LN.begin()); it != LN.end(); ++it ) if ( it->second > 0 ) isEnd = false;
		if  ( isEnd   ) { break; }
		if  ( toolong == 0 ) {
			cerr << "Error may happen in your bamfile " << bamInfile << endl;
			cerr << toolong << " lines has been lapsed, but I still cannot finish calculating the insertsize of RG: " << endl;
			for ( map<string, int>::iterator it(LN.begin()); it != LN.end(); ++it ) {
				if ( it->second > 0 ) cerr << it->first << endl;
			}
			cerr << "It's too long, may there something wrong! I have to stop the program and you have to chech your bamfile!" << endl;
			exit(1);
		}

		string rg, sa;
		if (!al.GetTag("RG", rg)) { 
            cerr << "Read " << al.Name << " doesn't have an alignment score RG:...  Cannot continue " << rg << endl;
        	exit(1);
        }
		if ( !al.IsProperPair()
			 || al.GetTag("SA", sa) 
			 || al.MapQuality < 30 // min mapped quality
			 || !al.IsMapped() 
			 || !al.IsMateMapped() 
			 || al.InsertSize < 0 ) continue;

		if ( !LN.count(rg) ) {
			cerr << "RG:Z:" << rg << " in alignment read " << al.Name << " does not appear on the header of your bamfile : " << bamInfile << endl;
		}

		ins[rg].push_back( al.InsertSize );
		--LN[rg];
		--toolong;
	}
	bamReader.Close();

	cerr << "# Read Group insertsize: Mean and SD\n";
	for ( map<string,vector<int> >::iterator it(ins.begin()); it != ins.end(); ++it ) { 
		is[it->first] = MeanAndStddev( it->second );
		cerr << "# " << it->first << "\t" << is[it->first].first << "\t" << is[it->first].second << "\n";
	}

	cerr << "\n***************** Calculaed the insert size Done *********************\n" << endl;
	return flag;
}

pair<int, int> MeanAndStddev ( vector<int> & v ) {

	pair<int, int> ms; // first->mean and second->sd
	double sum(0.0);
	// first work out the mean
	for ( size_t i(0); i < v.size(); ++i ) sum += v[i];
	ms.first = ceil( sum/(1.0 * v.size()) );
	
	sum = 0.0;
	for ( size_t i(0); i < v.size(); ++i ) sum += (1.0 * v[i] - ms.first) * (1.0 * v[i] - ms.first);
	ms.second = ceil( sqrt( sum / (1.0 * v.size()) ) );

	return ms;
}



