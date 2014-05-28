/* 
 * Author : Shujia Huang
 * Date   : 2013-10-2
 *
 * This is the header of API use for parse the *.axt file. 
 *
 */

#ifndef AXTVAR_H
#define AXTVAR_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include "Region.h"
#include "Fa.h"

using namespace std;

class Axt {
/*
1 chrM 16308 16389 Contig102837 1 81 + 5938
CATAGTACATAAAGTCATTTACCGTACATAGCACATTACAGTCAAATCCCTTCTCGTCCCCATGGATGACCCCCCTCAGATA
CATAGCACATATAGTCATTCATCGTACATAGCACATTATAGTCAAATCATTTCTCGTCCCCACGGAT-ATCCCCCTCAGATA

*/
public:
	string id;     // 1
    Region target; // chrM 16308 16389
    Region query;  // Contig102837 1 81
    char strand;   // +
    long score;    // 5938
    string tarSeq; // CATAGTACATAAAGTCATTTACCGTACATAGCACATTACAG (target sequence)
    string qrySeq; // CATAGCACATATAGTCATTCATCGTACATAGCACATTATAG (query  sequence)

public :
    Fa tarfa, qryfa;

public:

	void ConvQryCoordinate () {
	// This funtion just conversion the coordinate but not conversion the query sequence, 
	// which would be a problem when consider the sequence. So make this memerber function be private!
		if ( strand != '-' ) return;
		if ( qryfa.fa.empty() ) { 
			cerr << "[ERROR] Query fa sequence is missing! Make sure you've input it!" << endl; exit(EXIT_FAILURE);
		}

		if ( !qryfa.fa.count( query.id ) ) { 
			err ( "Missing some query id or query id can't match!!!\nThe unmatch query(main): " + query.id );
		}
		unsigned int itemp = query.start;
		query.start = qryfa.fa[query.id].length() - query.end + 1;
		query.end   = qryfa.fa[query.id].length() - itemp + 1;
	}

	void err(string reason) {
        cerr<<"\nError on ID: "<< id <<" ("<< target.id <<" <=> "<< query.id <<")\n";
        cerr<<"Reason => "     << reason   << endl;
		if ( SeqLength() > 100 ) { // Keep the output sequence is not too long!
			tarSeq = tarSeq.substr( 0, 100 ) + " ... "; 
			qrySeq = qrySeq.substr( 0, 100 ) + " ... ";
		}
		OutErrAlg (); exit(EXIT_FAILURE);
    }

	void CheckAxt() {
		RmGap2Gap();
		if ( tarSeq.length() != qrySeq.length() ) err ( "tarSeq.length() != qrySeq.length()" );
		if ( tarSeq.empty()  || qrySeq.empty()  ) err ( "tarSeq.empty()  || qrySeq.empty()"  );
		if ( strand != '+' && strand != '-'     ) err ( "strand != '+' && strand != '-'"     );
		if ( target.id.empty() || query.id.empty() ) err ( "target.id.empty() || query.id.empty()" );
		if ( target.start <= 0 || target.end <= 0 || query.start <= 0 || query.end <= 0 ) 
			err ( "target.start <= 0 || target.end <= 0 || query.start <= 0 || query.end <= 0" );
		if ( !qryfa.fa.empty() && !qryfa.fa.count( query.id ) )
			err ( "Missing some query id or query id can't match!!!\nThe unmatch query(main): " + query.id );
		if ( !tarfa.fa.empty() && !tarfa.fa.count( target.id ) )
			err ( "Missing some query id or target id can't match!!!\nThe unmatch query(main): " + target.id );
	}

	void OutErrAlg () { // Output the axt alignment to STDERR
		cerr<< "#"   << id << " " << target.id   << " " << target.start << " " << target.end << " "
        	<< query.id    << " " << query.start << " " << query.end    << " " << strand     << " " << score  << "\n"
         	<< tarSeq      << "\n"<< qrySeq      << "\n"<< endl;
	}

	void RmGap2Gap () { // Remove the gap to gap base into the algnment. it may calls by lastz's bug!
	// debug carefully here!!
		for ( string::iterator itt( tarSeq.begin() ), itq( qrySeq.begin() ); itt != tarSeq.end(); ) {
			if ( *itt == '-' && *itq == '-' ) {
				itt = tarSeq.erase( itt );
				itq = qrySeq.erase( itq );
			} else {
				++itt; ++itq;
			}
		}
	}
	unsigned int SeqLength () { return tarSeq.length(); }
};

// Call indel from per axt alignment.
#endif 








