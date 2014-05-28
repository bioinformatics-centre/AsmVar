#include "Stats.h"

short GetPairOrientation(BamAlignment& al, bool reverse) {

	short type;
    if (!(al.IsPaired() && al.IsMapped() && al.IsMateMapped())) {
        type =  UNPAIRED;
    } else if (al.RefID != al.MateRefID) {
       type =  DIFF_CHROM;
    } else if (al.IsReverseStrand() == al.IsMateReverseStrand()) {
        type = SAME;
    } else if ((al.Position <= al.MatePosition) == al.IsMateReverseStrand()) {
        type =  reverse ? OUTTER : INNER;
    } else if ((al.Position > al.MatePosition) == al.IsMateReverseStrand()) {
        type =  reverse ? INNER : OUTTER;
    } else {
        cerr << "Unexpected error in getPairOrientation().  Aborting." << endl;
        exit(1);
    }

	return type;
}

map< char, vector<pair<int32_t,int32_t> > > IndelSignalRefPos( int32_t pos, vector<CigarOp> cigarData ) {

	map< char, vector<pair<int32_t,int32_t> > > indel;
	int32_t mapLen = 0;
	for ( size_t i(0); i < cigarData.size(); ++i ) {
		switch ( cigarData[i].Type ) { // Refer: /path/bamtools/src/api/BamAlignment.cpp 'GetEndPosition' function
			case 'M' :
			case 'X' : 
			case 'N' :
			case 'P' : 
			case '=' : 
				mapLen += cigarData[i].Length ; break; 
			// Deletion signal cross reads
			case 'D' : 
				//indel['D'].push_back( make_pair(pos + mapLen - 1, pos + mapLen + cigarData[i].Length - 1) ); // Start and End
				indel['D'].push_back( make_pair(pos + mapLen, pos + mapLen + cigarData[i].Length - 1) ); // Start and End
				mapLen += cigarData[i].Length;
				break;
			case 'I' :
				indel['I'].push_back( make_pair(pos + mapLen - 1, pos + mapLen -1) );                        // Start and End

			// all other CIGAR chars do not affect end position
			default :
                break;
		}
	}
	return indel;
}

void LoadGC2Cov(string& filename, vector<float>& v_in) {

    ifstream ifs(filename.c_str());
    string line;

    if (!ifs.good()) {
        cerr << ERROR_PREFIX << "Error opening file '" << filename << "'" << endl;
        exit(1);
    }

    while (getline(ifs, line)) {
        vector<string> tmp;
        split("\t", line, tmp); //split

        size_t gc = atoi(tmp[0].c_str());

		// sanity check we've got the right GC
        if (gc != v_in.size()) {
            cerr << ERROR_PREFIX << "Error in GC to coverage file '" << filename << "'." << endl
                 << "Need GC in numerical order from 0 to 100.  Problem around this line:" << endl
                 << line << endl;
            exit(1);
        }

        v_in.push_back(atof(tmp[1].c_str()));
    }

    ifs.close();
}

void LoadGC(Tabix& ti, string& refname, vector<int>& v_in) {

    v_in.clear();          // The vector should be cleared!
    ti.setRegion(refname);
    vector<string> tmp;
    string line;

	// load the GC into vector and the position will be the index of vector +1 
	cerr << "** Load Reference "<< refname << " GC Contains" << endl;
    while (ti.getNextLine(line)) {
        split("\t", line, tmp);
        v_in.push_back(atoi(tmp[2].c_str()));
    }
}

void LoadGaps(string fname, map<string, list<pair<long, long> > >& gaps) {

    Tabix ti(fname);
    string line;
	
	cerr << "\n*****************  Loading reference gaps' regions  *****************" << endl;
    while (ti.getNextLine(line)) {
        vector<string> data;
        split("\t", line, data);
        gaps[data[0]].push_back( make_pair(atoi(data[1].c_str()), atoi(data[2].c_str()) ) );
    }
	cerr << "*****************        Reference gaps loaded      *****************\n" << endl;
}

bool CheckStatsAndFitLength ( Stats& stats ) { // Just be used in the OutputStats function

	bool isGood(true);
	vector<long> positions;
	positions.push_back( stats.properPair.Position()       );
	positions.push_back( stats.pLowAlignScore.Position()   );
	positions.push_back( stats.crossSignal.Position()      );
	positions.push_back( stats.badInserSizePair.Position() );
	positions.push_back( stats.wrongOrientation.Position() );
	positions.push_back( stats.clipAlign.Position()        );
	positions.push_back( stats.single.Position()           );
	
	if ( positions.empty() ) { 
		isGood = false;
	} else {
		long p = positions[0];
		for ( size_t i(1); i < positions.size(); ++i ) {
			if ( p != positions[i] ) isGood = false;
		}
	}
	return isGood;
}

void OutputStats (string refname, long pos, Stats& stats, map<string, list<pair<long, long> > >& gaps) {

	// see if any gaps in this reference sequence
	map<string, list<pair<long, long> > >::iterator gapsIter = gaps.find(refname);
	list<pair<long, long> >::iterator gIter;
	// if do have gaps!
	if ( gapsIter != gaps.end() ) gIter = gapsIter->second.begin();

	if ( !CheckStatsAndFitLength( stats ) ) {
		cerr << ERROR_PREFIX << "the values 'stats' looks not good. They haven't the same coordinates in coverage calculate.\n";
		exit(1);
	}
	long p = stats.properPair.Position();
	// Output all the covered stats of each position left to 'pos' not include 'pos'
	for ( ; p < pos; ++p ) {

		bool inGap(false);
		if ( gapsIter != gaps.end() ) {
			while ( gIter != gapsIter->second.end() && gIter->second < p ) ++gIter;
			if    ( gIter != gapsIter->second.end() && gIter->first <= p ) inGap = true;
		}
		//float expectedCovDep = (inGap) ? 0 : stats.gc2cov[ stats.refGC[p-1] ];
		//float correctFactor  = (inGap) ? 0 : 1.0 *   calculate this factor later.
		
		long allcov = stats.properPair.Depth() + stats.pLowAlignScore.Depth() + stats.crossSignal.Depth() + stats.badInserSizePair.Depth() + 
					  stats.wrongOrientation.Depth() + stats.single.Depth()   + stats.clipAlign.Depth();

		if ( p != stats.properPair.Position() ) cerr << ERROR_PREFIX << "Position Error. " << stats.properPair.Position() << " != " << p << endl;
		/*
		cout << refname << "\t" << p << "\t" << allcov                           << "\t"
			 << (allcov > 0 ? 1.0 * stats.properPair.Depth()       / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.pLowAlignScore.Depth()   / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.crossSignal.Depth()      / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.badInserSizePair.Depth() / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.wrongOrientation.Depth() / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.single.Depth()           / allcov : 0) << "\t"
			 << (allcov > 0 ? 1.0 * stats.clipAlign.Depth()      / allcov : 0) << "\t-" << endl;
			 //<< stats.clipAlign.Depth()  << "\t" << expectedCovDep             << endl;
		*/
		cout << refname << "\t" << p << "\t" << allcov << "\t"
			 << stats.properPair.Depth()               << "\t"
			 << stats.pLowAlignScore.Depth()           << "\t"
			 << stats.crossSignal.Depth()              << "\t"
			 << stats.badInserSizePair.Depth()         << "\t"
			 << stats.wrongOrientation.Depth()         << "\t"
			 << stats.single.Depth()                   << "\t"
			 << stats.clipAlign.Depth() << endl;
			 //<< stats.clipAlign.Depth()  << "\t" << expectedCovDep             << endl;

		// Move to the next position
		stats.properPair.Next();
    	stats.pLowAlignScore.Next();
		stats.crossSignal.Next();
		stats.badInserSizePair.Next();
		stats.wrongOrientation.Next();
		stats.clipAlign.Next();
    	stats.single.Next();
	}
}







