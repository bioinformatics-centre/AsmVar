// base.inl
// Ver 1.0
// Author : Shujia Huang
// Date   : 2011/04/06

#include "Join.h"

string join ( const char *delim, vector<string>& goven ) {

    string strLine;
    strLine.clear();

    if ( goven.empty() ) return strLine;

    strLine = goven[0];
    for ( size_t i(1); i < goven.size(); ++i ) {

        strLine += (delim + goven[i] );
    }

    return strLine;
}

string join ( char* goven[], int size ) {

	vector< string > str_goven;
	for (size_t i (0); i < size; ++i ) {
		str_goven.push_back( goven[i] );
	}

	return join ( " ", str_goven );
}

//2011/08/23
string toupper ( string strLine ){ 

	transform( strLine.begin(), strLine.end(), strLine.begin(), (int (*)(int) ) toupper ); 
	return strLine;
}

// Calculate localtime
char* local_time () {

	time_t rawtime;
	time ( &rawtime );
	return ( asctime(localtime(&rawtime)) );
}




