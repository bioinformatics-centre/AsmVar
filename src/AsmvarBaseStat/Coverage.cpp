#include "Coverage.h"

Coverage::Coverage() {

	Coverage( 0, 0 );
}

Coverage::Coverage(long n, long pos) : coord ( pos ) {

	covdep = deque<long>( n+1, 0 ); // Add 1 to avoid empty
}

Coverage::Coverage(const Coverage &C ) {

	coord = C.coord;
	covdep= C.covdep;
}

void Coverage::Next () {

	++coord;
	covdep.pop_front( );
	covdep.push_back(0);
}

long Coverage::Position() {
	return coord;
}

long Coverage::Depth () {

	return covdep[0];
}

size_t Coverage::Size() {

	return covdep.size();
}

void Coverage::Add( long start, long end ) {

	//if ( start > end ) { cerr << "Error in Coverage::Add( long start, long end ). The start > end\n" << start << " > " << end; exit(1); }
	if ( start > end ) { 
		cerr << "Warning in Coverage::Add( long start, long end ). The start > end\n" << start << " > " << end;
		if ( start - end == 1 ) {
			end = start;
			cerr << " ... Continue ... " << endl;
		} else {
			cerr << " ... Cannot Continue ... " << endl; exit(1); 
		}
	}

	for ( long p(start); p < end + 1; ++p ) Add ( p );
}

void Coverage::Add( long pos ) {

	if (pos < coord) { 
		cerr <<"Coverage ERROR. Do not allow the position you input ( "      << pos 
             <<" ) is less than the start (" << coord << ") of this region." << endl;
		exit(1);
	}
	if (covdep.size() + coord < pos + 1) { // Shift the size of covdep when position overflow 
		long remain = pos - (covdep.size()+coord - 1);
		for ( long i(0); i < remain; ++i ) covdep.push_back(0);
	}

	++covdep[pos-coord];
}

