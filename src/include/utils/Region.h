/* Author : Shujia Huang
 * Date   : 2013-10-2
 *
 * Region class
 */
#ifndef __REGION_H__
#define __REGION_H__

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h> 

using namespace std;

class Region {

public:
    string id;
    long int start;
    long int end;
	long int regsize; // The region size
	string info;

public:
    Region () : start(0), end (0), regsize(0) { id.clear(); info.clear(); }
    Region (const Region &R) { // The copy construct function. 
		id = R.id; start = R.start; end = R.end; info = R.info; 
		regsize = R.regsize;
	}
	void SetRegSize() { regsize = (start + end > 0) ? labs(end-start) + 1 : 0; }

public:

	void OutErrReg() {
		cerr << "# " << id << "\t" << start << "\t" << end << "\t" << info << "\n";
    }

	bool IsOverlap (Region &R) { // determine this region overlap with other region or not!
		bool flag (false);
		if ( (id == R.id) && (start <= R.end && end >= R.start) ) flag = true;
		return flag;
	}
	bool isEmpty () { return (end == 0); }
};

vector<Region> MergeRegion(vector<Region> &RegVector, int delta);

#endif

