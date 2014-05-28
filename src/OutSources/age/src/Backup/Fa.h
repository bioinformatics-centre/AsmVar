/* Author : Shujia Huang
 * Date   : 2013-10-10
 * 
 */
#ifndef FA_H
#define FA_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdlib.h>
#include "gzstream.h"

using namespace std;
typedef unsigned int uint; // The biggest number of uint is 2^32 = 4294967296

class Fa {

public:

	map<string, string> fa;
	uint length; // Recorde the total size 0f fa. Should not bigger than 4G!
public:
    Fa() : length(0) { fa.clear(); }

    void Load( string file ) { 
        igzstream I( file.c_str() );
        if ( !I ){
            cerr << "Cannot open file : " << file.c_str() << endl;
            exit(1);
        }

        string tmp, faId;
        while ( 1 ) {

            I >> tmp;
            if ( I.eof() ) break;
            if ( tmp[0] != '>' ) {

                fa[faId] += tmp;
				length   += tmp.length();
            } else {
                faId.assign ( tmp, 1, string::npos );
                fa[faId].clear();
            }
            getline ( I, tmp, '\n' );
        }
        I.close();
    }

	uint Nlength( string id, uint start, uint end ) {
		if ( end < start )    { cerr << "[ERROR] start > end, in Nlength( string id, uint start, uint end ) in class 'Fa'\n"; exit(1); }
		if ( !fa.count( id ) ){ cerr << "[ERROR] There is no " << id << " in the fasta reference.\n" ; exit(1); }
		string str = fa[id].substr( start - 1, end - start + 1 );

		uint len(0);
		for ( size_t i(0); i < str.size(); ++i ) {
			if ( str[i] == 'N' || str[i] == 'n' ) ++len;
		}
		return len;
	}

	void Code () {
		// assign the characters to 1, instead of 0, cause if 0,it'll be ambiguous in (find_first_not_of(0) in Output() function.).
		for ( map<string, string>::iterator it( fa.begin() ); it != fa.end(); ++it ) it->second.assign( it->second.length(), 1 ); 
	}
};

#endif

