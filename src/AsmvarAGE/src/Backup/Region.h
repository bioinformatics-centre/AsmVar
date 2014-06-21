/* Author : Shujia Huang
 * Date   : 2013-10-2
 *
 * Region class
 */
#ifndef REGION_H
#define REGION_H

#include <iostream>

using namespace std;

class Region {

public:
    string id;
    unsigned long int start;
    unsigned long int end;
	string info;

public:
    Region () : start(0), end (0) { id.clear(); info.clear(); }
    Region ( const Region & R ) { id = R.id; start = R.start; end = R.end; info = R.info; } // The copy construct function.

public:
	bool IsOverlap ( Region & R ) { // determine overlap other region or not!
		bool flag (false);
		if ( (id == R.id) && (start <= R.end && end >= R.start) ) flag = true;
		return flag;
	}
	bool isEmpty () { return (end == 0); }
};

#endif

