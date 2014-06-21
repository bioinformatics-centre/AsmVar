// base.h
// Ver 1.0
// Author : Shujia Huang
// Date   : 2011/04/06

#ifndef JOIN_H
#define JOIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>

using namespace std;

string join ( const char *delim, vector<string>& goven );
string join ( char* goven[], int size ); // New 2013-05-28 17:20:18
string toupper ( string strLine );
//Calculation locaktime
char* local_time ();

#endif

