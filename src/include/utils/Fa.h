/* Author : Shujia Huang
 * Date   : 2013-10-10
 * 
 */
#ifndef __FA_H__
#define __FA_H__

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdlib.h>
#include "gzstream.h"

using namespace std;

//typedef unsigned long int uint; // The biggest number of uint is 2^32 = 4294967296

class Fa {

public:

	map<string, string> fa;
	map<int, string> order2id; // Record the order of fa id
	uint length; // Recorde the total size of fa. Should not bigger than 4G!
	uint nsize;  // Recorde the total size of 'N' in fa. Should not bigger than 4G!

public:
    Fa() : length(0), nsize(0) { fa.clear(); }
	Fa( string file ) { Load(file); }

    void Load( string file ) { 

        igzstream I( file.c_str() );
        if ( !I ){
            cerr << "Cannot open file : " << file.c_str() << endl;
            exit(1);
        }

        string tmp, faId;
		int order = 1;
        while ( 1 ) {

            I >> tmp;
            if ( I.eof() ) break;
            if ( tmp[0] != '>' ) {

                fa[faId] += tmp;
				length   += tmp.length();
				for (size_t i(0); i < tmp.size(); ++i) 
					if (tmp[i] == 'N' || tmp[i] == 'n') ++nsize;

            } else {
                faId.assign ( tmp, 1, string::npos );
                if (fa.count(faId)) {
                    cerr << "[ERROR] Multiple fa id(" << faId 
                         << ") in file: " << file << "\n" << endl;
                    exit(1);
                }
                fa[faId].clear();

				order2id[order] = faId;
				++order;
            }
            getline ( I, tmp, '\n' );
        }
        I.close();
		cerr << "[INFO] " << file << " complete loading." << endl;
    }

public :
	void Clear () { fa.clear(); length = 0; }

	void CheckFaId(string id) {

		if (!fa.count(id)) {
			cerr << "[ERROR] Missing some fa id or fa id can't match in "
				 << "\nThe unmatch id: " + id << "\n";
			exit(1);
		}
	}

	int Nlength(string id, uint start, uint end) {
		if ( end < start )    { cerr << "[ERROR] start > end" << start << ", " << end << " in Nlength( string id, uint start, uint end ) in class 'Fa'\n"; exit(1); }
		if ( !fa.count( id ) ){ cerr << "[ERROR] There is no " << id << " in the fasta reference.\n" ; exit(1); }
		if ( end > fa[id].length() ) end = fa[id].length();
		string str = fa[id].substr( start >= 1 ? start - 1: 0, end - start + 1 );

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

