#ifndef COVERAGE_H
#define COVERAGE_H

#include <iostream>
#include <deque>
#include <stdlib.h>     /* exit, EXIT_FAILURE */

using namespace std;

class Coverage {

public:
	Coverage();
	Coverage(long n, long pos);
	Coverage(const Coverage &C); // Copy constructor function

	// returns the position of current base
	long Position();

	// returns the coverage at the current base
	long Depth();

	// increase the depth from position 'start' to position 'end'
	void Add( long start, long end );

	// add a read that ends at position pos
	void Add (long pos);

	// move the coord 1 base to the right
	void Next ();

	// return the size of 'covdep'
	size_t Size();
private:

	long coord;         // Just mark the coordinate covdep[0].
	deque<long> covdep; // The first position of this deque is 'coord'
};

#endif

