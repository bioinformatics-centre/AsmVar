// base.inl
// Ver 1.0
// Author : Shujia Huang
// Date   : 2011/04/06

#include "utility.h"

// split Function
void split( const char *delim, string strLine, vector<string>& tokens, bool is_add ) {

	if ( !is_add ) {
		tokens.clear();
	}
	while( 1 ) {
		//erase delimiter
	   int i = strLine.find_first_not_of(delim);

	   if(i == -1) break;

	   strLine.erase(0, i);

	   i = strLine.find_first_of(delim);
	   if(i == -1) {

		    tokens.push_back(strLine);
		    break;
	   } else {
		    string token = strLine.substr(0, i);
		    strLine.erase(0, i);
		    tokens.push_back(token);
	   }
	}
}

// 2011/08/23
void split( const char *delim, string strLine, vector<char>& tokens, bool is_add ) {

	if ( !is_add ) tokens.clear();
	while ( 1 ) {

		//erase delimiter
		int i = strLine.find_first_not_of(delim);
		if ( i == -1 ) break;

		strLine.erase(0, i);

		i = strLine.find_first_of(delim);
		if ( i == -1 ) {
			tokens.push_back ( strLine[0] ); // string -> char
			break;
		} else {
			string token = strLine.substr(0, i);
			if ( token.length() > 1 ) {
				cerr << "ERROR:String element cannot put in char vector in split function." << endl;
				exit(1);
			}
			strLine.erase(0, i);
			tokens.push_back( token[0] );
		}

	}
}
// 2011/08/23
void split( const char *delim, string strLine, vector<int>& tokens, bool is_add ) {

	if ( !is_add ) tokens.clear();

	vector<string> str_int;
	split ( delim, strLine, str_int, is_add );

	for ( size_t i(0); i < str_int.size(); ++i ) {
		tokens.push_back( atoi( str_int[i].c_str() ) );
	}
}

void split( const char *delim, string strLine, list<string>& tokens, bool is_add ) {

    if ( !is_add ) {
        tokens.clear();
    }
    while( 1 ) {
        //erase delimiter
       int i = strLine.find_first_not_of(delim);

       if(i == -1) break;

       strLine.erase(0, i);

       i = strLine.find_first_of(delim);
       if(i == -1) {

            tokens.push_back(strLine);
            break;
       } else {
            string token = strLine.substr(0, i);
            strLine.erase(0, i);
            tokens.push_back(token);
       }
    }
}

void split( const char *delim, string strLine, list<char>& tokens, bool is_add ) {

    if ( !is_add ) tokens.clear();
    while ( 1 ) {

        //erase delimiter
        int i = strLine.find_first_not_of(delim);
        if ( i == -1 ) break;

        strLine.erase(0, i);

        i = strLine.find_first_of(delim);
        if ( i == -1 ) {
            tokens.push_back ( strLine[0] ); // string -> char
            break;
        } else {
            string token = strLine.substr(0, i);
            if ( token.length() > 1 ) {
                cerr << "ERROR:String element cannot put in char vector in split function." << endl;
                exit(1);
            }
            strLine.erase(0, i);
            tokens.push_back( token[0] );
        }

    }
}

void split( const char *delim, string strLine, list<int>& tokens, bool is_add ) {

    if ( !is_add ) tokens.clear();

    vector<string> str_int;
    split ( delim, strLine, str_int, is_add );

    for ( size_t i(0); i < str_int.size(); ++i ) {
        tokens.push_back( atoi( str_int[i].c_str() ) );
    }
}


void join ( const char *delim, vector<string>& goven, string &strLine ) {

	strLine.clear();
	if ( goven.empty() ) return;

	strLine = goven[0];
	for ( size_t i(1); i < goven.size(); ++i ) {

		strLine += (delim + goven[i] );
	}
}


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

// 2011/08/23
string join ( const char *delim, vector<char>& goven ) {
	
	string strLine;
	strLine.clear();

	if ( goven.empty() ) return strLine;
	strLine.append( char2str(goven[0]) );
	for ( size_t i(1); i < goven.size(); ++i ) {

		strLine.append ( delim + char2str( goven[i] ) );
	}

	return strLine;
}
//2011-05-21
string join ( const char *delim, vector<int>& goven ) {

	vector<string> str_goven;
	for ( size_t i(0); i < goven.size(); ++i ) str_goven.push_back( itoa( goven[i] ) );

	return join( delim, str_goven );
}

string join ( const char *delim, list<string> &goven ){

	string strLine;
	strLine.clear();

	if ( goven.empty() ) return strLine;

	strLine = goven.front();
	list<string>::iterator iter = goven.begin();
	for ( ++iter; iter != goven.end(); ++iter ) strLine += ( delim + *iter );

	return strLine;
}

string join ( const char *delim, list<char> &goven ){

	string strLine;
    strLine.clear();

    if ( goven.empty() ) return strLine;
	list< char >::iterator iter = goven.begin();
    strLine = char2str( *iter );
	for ( ++iter; iter != goven.end(); ++iter ) strLine += ( delim + char2str( *iter ) );

    return strLine;
}

string join ( const char *delim, list<int> &goven ){

	vector<string> str_goven;
	list<int>::iterator iter = goven.begin();
	for (; iter != goven.end(); ++iter ) str_goven.push_back( itoa( *iter ) );

	return join( delim, str_goven );
}

string join ( const char *delim, deque<string> &goven ) {

	string strLine;
    strLine.clear();

    if ( goven.empty() ) return strLine;

    strLine = goven[0];
    for ( size_t i(1); i < goven.size(); ++i ) {

        strLine += (delim + goven[i] );
    }

    return strLine;
}

string join ( const char *delim, deque<char> &goven ) {

	string strLine;
    strLine.clear();

    if ( goven.empty() ) return strLine;
    strLine.append( char2str(goven[0]) );
    for ( size_t i(1); i < goven.size(); ++i ) {

        strLine.append ( delim + char2str( goven[i] ) );
    }

    return strLine;
}

string join ( const char *delim, deque<int> &goven) {

	deque<string> str_goven;
    for ( size_t i(0); i < goven.size(); ++i ) str_goven.push_back( itoa( goven[i] ) );

    return join( delim, str_goven );
}

string join ( char* goven[], int size ) {

	vector< string > str_goven;
	for (size_t i (0); i < size; ++i ) {
		str_goven.push_back( goven[i] );
	}

	return join ( " ", str_goven );
}

//2011-05-21
/*
string itoa ( const int number ) {

	stringstream num2str;
	num2str << number;

	return num2str.str();
}
*/
string itoa( int number ) {

    vector< int > meta;
	bool isnegtive = (number < 0); // Negative number. e.g : -22
	number = (number < 0) ? -number : number;
    do {
        meta.push_back( number % 10 );
    } while ( number /= 10 );

    string str = (isnegtive) ? "-" : "";
    while ( !meta.empty() ) {

        str.insert( str.end(), meta.back() + '0' );
        meta.pop_back();
    }

	return str;
}

string ftoa ( double dbl ) {

	ostringstream strs;
	strs << dbl;
	return strs.str();
}


//2011/08/23
string char2str( char c ) {
	
	string str;
	str.insert( str.begin(), c );

	return str;
}
string toupper ( string strLine ){ 

	transform( strLine.begin(), strLine.end(), strLine.begin(), (int (*)(int) ) toupper ); 
	return strLine;
}

//
char mytoupper( const char c ) {

    if ( c < 'A' || c > 'z' || ( c > 'Z' && c < 'a' ) ) return c;
    return c & 0x5f;
}

char mytolower( const char c ) {
    if ( c < 'A' || c > 'z' || ( c > 'Z' && c < 'a' ) ) return c;
    return c | 0x20;
}


void open_file ( ifstream &in, const string &file_name ) {

	in.close();
	in.clear();
	
	in.open ( file_name.c_str() );
}

// Calculate localtime
char* local_time () {

	time_t rawtime;
	time ( &rawtime );
	return ( asctime(localtime(&rawtime)) );
}


// Math function
template < class T >
double mean( vector<T> &values ) {

	if ( values.empty() ) {

		cerr << "ERROR:The vector is empty. It cannot use to calculate variance" << endl;
		exit(1);
	}

	double total = 0;
	for ( size_t i(0); i < values.size(); ++i ) {
		total += values[i];
	}

	return total / values.size();
}

template < class T >
double mean( deque<T> &values ) {

	if ( values.empty() ) {

		cerr << "ERROR:The deque is empty. It cannot use to calculate mean." << endl;
		exit(1);
	}

	double total = 0;
	for ( size_t i(0); i < values.size(); ++i ) {
		total += values[i];
	}

	return total / values.size();
}

template < class T >
double variance ( vector<T> &values ) {

	double mean_value = mean ( values );
	double var = 0;

	if ( values.empty() ) {

		cerr << "ERROR:The vector is empty. It cannot use to calculate variance." << endl;
		exit(1);
	}

	for ( size_t i(0); i < values.size(); ++i ) var += pow ( (values[i] - mean), 2 );

	return var / values.size();
}

template < class T >
double variance ( deque<T> &values ) {

	double mean_value = mean ( values );
	double var = 0;

	if ( values.empty() ) {

		cerr << "ERROR:The deque is empty. It cannot use to calculate variance" << endl;
		exit(1);
	}

	for ( size_t i(0); i < values.size(); ++i ) var += pow ( (values[i] - mean_value ), 2 );

	return var /  values.size();
}

//Trim spaces among string  2011/09/28
string trim_space( string str ) {

	string::iterator it = str.begin();
	for ( ; it != str.end(); ) {

		if ( (*it) == ' ' || (*it) == '\t' || (*it) == '\n' ) { 
			str.erase( it );
		} else {
			++it;
		}
	}

	return str;
}





