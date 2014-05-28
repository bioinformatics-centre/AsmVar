/****************************************************************************
 *
 *     AGE -- Alignment with Gap Excision
 *     Copyright (C)  Alexej Abyzov
 *                
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Creative Commons license
 * (Attribution-NonCommerical).
 * See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
 *
 * Author: Alexej Abyzov
 */

#ifndef __SEQUENCE_HH__
#define __SEQUENCE_HH__

// Nucleotide sequence
// Allowed nucleotides are A,C,T,G,U,X,N in both lower and upper cases
// Other nucleotides (R,Y,M,K,S,W,B,D,H,V) are also alowed by not scored
class Sequence
{
private:
	string _seq,_name;
	int    _start;     // 1 based coordinate
	bool   _reverse;

public:
	inline string &name()      { return _name; }
	inline string &sequence()  { return _seq; }
	inline int     start()     { return _start; }  // 1 based coordinate
	inline int     reverse()   { return _reverse; }

public:
	Sequence(const char *seq) : _seq(seq), _name(""),_start(1),_reverse(false), _next(NULL),_prev(NULL) {
		string ret = "";
		for (int i = 0;i < _seq.length();i++) {
			char c = _seq[i];
			if (c == ' ') continue;      // Space
			else if (isNuc(c)) ret += c; // Nucleotide
			else if (!isGap(c)) {        // Gap
				cerr<<"Unrecognized nucleotide '"<<c<<"' in sequence."<<endl;
				_seq = "";
				return;
			}
		}
		_seq = ret;
	}

	Sequence(string file,int start = 1,int end = -1) : _seq(""), _name(""), _start(1),_reverse(false), _next(NULL),_prev(NULL) {

		if (start > 0 && end > 0 && end < start) {
			cerr<<"Invalid coordinates ("<<start<<", "<<end<<"). "
				<<"Start is larger than end."<<endl;
			cerr<<"No sequence read!"<<endl;
			return;
		}

		if (start >= 0) _start = start;

		ifstream in(file.c_str());
		if (!in) {
			cerr<<"Can't open input file '"<<file<<"'."<<endl;
			return;
		}

		char l[1024];
		while (!in.eof() && in.getline(l,1024)) 
			if (in.gcount() != 0) break;

		if (l[0] != '>') {
			cerr<<"Can't find fasta header in the file '"<<file<<"'."<<endl;
			return;
		}
		_name = &l[1];

		string line = "";
		int n_passed = 0;
		while (!in.eof()) {
			in>>line;
			string checked = "";
			for (int i = 0;i < line.length();i++) {
				char c = line[i];
				if (c == ' ') continue;          // Space
				else if (isNuc(c)) checked += c; // Nucleotide
				else if (!isGap(c)) {            // Gap
					cerr<<"Unrecognized nucleotide '"<<c
						<<"' in file '"<<file<<"'."<<endl;
					_seq = "";
					return;
				}
			}
			int n_checked = checked.length();
			if (n_checked == 0) break;

			int new_passed = n_passed + n_checked;
			if (new_passed < start) { ; // Not reached start yet
			} else if (n_passed < start && new_passed >= end) { // Spans whole region
				_seq += checked.substr(start - n_passed - 1,end - start + 1);
			} else if (n_passed < start && new_passed >= start) { // Spans start only
				_seq += checked.substr(start - n_passed - 1);
			} else if (end < start) { // No end
				_seq += checked;
			} else if (new_passed < end) { // Not reached end yet
				_seq += checked;
			} else if (n_passed < end && new_passed >= end) { // Spans end only
				_seq += checked.substr(0,end - n_passed);
			}
			n_passed = new_passed;
			line = "";
		}
		in.close();
	}

public :
	Sequence(string &seq, string &name, int start, bool rev) : _seq(seq),
		_name(name),
		_start(start),
		_reverse(rev),
		_next(NULL),
		_prev(NULL)
	{}
	Sequence ( const Sequence& S ) {
		this->_seq     = S._seq;
		this->_name    = S._name;
		this->_reverse = S._reverse;
		this->_next    = S._next;
		this->_prev    = S._prev;
	}
	~Sequence () {
	// I should add code here! And not just use 'deleteSequences()' to release my object
	}

public :
	Sequence *clone() {
		return new Sequence(_seq,_name,_start,_reverse);
	}

	static char complement(char c) {
		if      (c == 'a') return 't';
		else if (c == 'c') return 'g';
		else if (c == 't') return 'a';
		else if (c == 'g') return 'c';
		else if (c == 'A') return 'T';
		else if (c == 'C') return 'G';
		else if (c == 'T') return 'A';
		else if (c == 'G') return 'C';
		else if (c == 'u') return 'a';
		else if (c == 'U') return 'A';
		return c;
	}

	static bool isNuc(char c) {
		if (c == 'A' || c == 'a' ||
			c == 'T' || c == 't' ||
			c == 'C' || c == 'c' ||
			c == 'G' || c == 'g' ||
			c == 'U' || c == 'u' ||
			c == 'X' || c == 'x' ||
			c == 'N' || c == 'n' ||
			c == 'R' || c == 'r' ||
			c == 'Y' || c == 'y' ||
			c == 'M' || c == 'm' ||
			c == 'K' || c == 'k' ||
			c == 'S' || c == 's' ||
			c == 'W' || c == 'w' ||
			c == 'B' || c == 'b' ||
			c == 'D' || c == 'd' ||
			c == 'H' || c == 'h' ||
			c == 'V' || c == 'v') return true;
		return false;
	}

	static char gap() { return '-'; }

	static bool isGap(char c) {
		//if (c == '-') return true;
		//if (c == ' ') return true;
		//return false;
		return (c == ' ' || c == '-'); // By Shujia Huang 2014-04-08 20:27:55
	}

	static bool sameNuc(char a,char b) {
		char aa = toupper(a), bb = toupper(b);
		if ((aa == 'A' || aa == 'C' || aa == 'T' || aa == 'G' || aa == 'U') &&
			(bb == 'A' || bb == 'C' || bb == 'T' || bb == 'G' || bb == 'U') &&
			aa == bb) return true;
		return false;
	}

	void revcom() {

		if (_reverse) _start -= _seq.length() - 1;
		else          _start += _seq.length() - 1;
		_reverse = !_reverse;
		string ret = "";
		for (int i = _seq.length() - 1;i >= 0;i--) ret += complement(_seq.at(i));
		_seq = ret;
	}

	Sequence *substr(int start,int end) { // 1 based coordinates
	
		if (start > 0 && end > 0 && start > end) return NULL;
		if (start <= 0)                        start = 1;
		if (end   <= 0 || end > _seq.length()) end   = _seq.length();
		int s = start, e = end;
		if (_reverse) {
			//s = _seq.length() - e + 1;
			//e = _seq.length() - s + 1; // This is not Right!!!! By Shujia Huang 2014-04-07 15:32:33
			s = _seq.length() - end   + 1; // By Shujia Huang 2014-04-07 15:32:33
			e = _seq.length() - start + 1; // By Shujia Huang 2014-04-07 15:32:33
		}
		string tmp = "";
		if (_seq.length() > 0) tmp = _seq.substr(s - 1,e - s + 1);
		return new Sequence(tmp,_name,_start + start - 1,_reverse);
	}

	// Next and previous sequences
private:
	Sequence *_next,*_prev;

public:
	Sequence *next() { return _next; }
	Sequence *prev() { return _prev; }

	// Adding sequence before
	bool addBefore(Sequence *newSequence) {
		// Check if input is good                                                
		if (!newSequence) return false;

		// Check if same object                                                  
		if (newSequence == this) return false;

		// Check if the call comes from addAfter                                 
		if (newSequence->next() == this && !_prev) {
			_prev = newSequence;
			return true;
		}

		// Chack if they are already paired                                      
		if (newSequence->prev() == this && _next == newSequence) return true;

		// Check if it can be added                                              
		if (_prev) return false;
		if (newSequence->next()) return false;
		if (newSequence->prev()) return false;

		_prev = newSequence; // Set previous

		if (!newSequence->addAfter(this)) { // Update next for newSequence
			_prev = NULL;
			return false;
		}
		return true;
	}

	// Adding atom after                                                         
	bool addAfter(Sequence *newSequence) {
		// Check if not null                                                     
		if (!newSequence) return false;

		// Check if same object                                                  
		if (newSequence == this) return true;

		// Check if the call comes from addAfter                                 
		if (newSequence->prev() == this && !_next) {
			_next = newSequence;
			return true;
		}

		// Chack if they are already paired                                      
		if (newSequence->prev() == this && _next == newSequence) return true;

		// Check if it can be added                                              
		if (_next) return false;
		if (newSequence->next()) return false;
		if (newSequence->prev()) return false;

		_next = newSequence; // Set next

		if (!newSequence->addBefore(this)) { // Update previous for newSequence
			_next = NULL;
			return false;
		}
		return true;
	}

	static Sequence* parseSequences(string file) {

		Sequence *first = NULL,*last = NULL;
		ifstream in(file.c_str());
		if (!in) {
			cerr<<"Can't open input file '"<<file<<"'."<<endl;
			return first;
		}

		string name = "",seq = "";
		char *line  = new char[65536];
		while (!in.eof()) {

			in.getline(line,65536);
			int n = in.gcount();
			if (n == 0) continue;

			if (line[0] == '>') {
				if (name.length() > 0) {
					Sequence *s = new Sequence(seq,name,1,false);
					if (first) {
						last->addAfter(s);
						last = s;
					} else first = last = s;
				}
				name = &line[1];
				seq  = "";
				continue;
			}

			string checked = "";
			for (int i = 0;i < n - 1;i++) {
				char c = line[i];
				if (isNuc(c)) checked += c;
				else if (!isGap(c)) {
					cerr<<"Unrecognized nucleotide '"<<c
						<<"' in file '"<<file<<"'."<<endl;
					seq  = "";
					name = "";
				}
			}
			seq += checked;
		}
		in.close();

		if (name.length() > 0) {
			Sequence *s = new Sequence(seq,name,1,false);
			if (first) {
				last->addAfter(s);
				last = s;
			} else first = last = s;
		}

		delete[] line;
		return first;
	}

	static void deleteSequences(Sequence *seqs) {
		Sequence *s = seqs;
		while (s) {
			Sequence *next = s->next();
			delete s;
			s = next;
		}
	}
};

#endif

