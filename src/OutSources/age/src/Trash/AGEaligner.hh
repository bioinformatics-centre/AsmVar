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

#ifndef __AGEALIGNER_HH__
#define __AGEALIGNER_HH__

//--- C/C++ includes ---
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#ifdef OMP
#include <omp.h>
#endif
#ifdef AGE_TIME
#include <sys/time.h>
#endif
using namespace std;

//--- AGE includes ---
#include "Sequence.hh"
#include "Scorer.hh"

const static unsigned char DIAGONAL   = 1;
const static unsigned char HORIZONTAL = 2;
const static unsigned char VERTICAL   = 3;
const static unsigned char MASK       = 3;

const static unsigned char FS_DIAGONAL   = DIAGONAL;
const static unsigned char FS_HORIZONTAL = HORIZONTAL;
const static unsigned char FS_VERTICAL   = VERTICAL;
const static unsigned char FS_MASK       = MASK;

const static unsigned char RS_DIAGONAL   = DIAGONAL<<2;
const static unsigned char RS_HORIZONTAL = HORIZONTAL<<2;
const static unsigned char RS_VERTICAL   = VERTICAL<<2;
const static unsigned char RS_MASK       = MASK<<2;

const static unsigned char FM_DIAGONAL   = DIAGONAL<<4;
const static unsigned char FM_HORIZONTAL = HORIZONTAL<<4;
const static unsigned char FM_VERTICAL   = VERTICAL<<4;
const static unsigned char FM_MASK       = MASK<<4;

const static unsigned char RM_DIAGONAL   = DIAGONAL<<6;
const static unsigned char RM_HORIZONTAL = HORIZONTAL<<6;
const static unsigned char RM_VERTICAL   = VERTICAL<<6;
const static unsigned char RM_MASK       = MASK<<6;

class Block;
class AliFragment;

class AGEaligner
{
public:
  const static int INDEL_FLAG        = 0x01;
  const static int TDUPLICATION_FLAG = 0x02;
  const static int INVERSION_FLAG    = 0x04;
  const static int INVR_FLAG         = 0x08;
  const static int INVL_FLAG         = 0x10;
private:
  const static int ALL_FLAGS         = INDEL_FLAG | TDUPLICATION_FLAG |
    INVERSION_FLAG | INVR_FLAG | INVL_FLAG;

private:
  Sequence       &_s1,  *_s1_rc,  &_s2;    // Sequences with headers
  string         &_seq1,&_seq1_rc,&_seq2;  // Nucleotide sequences
  int             _len1,           _len2;  // Lengths
  AliFragment    *_frags,*_frags_alt;      // Fragment of alignment
  unsigned short *f_score,*r_score;        // Scoring matrices
  unsigned char  *trace;                   // Trace back matrix
  int            score_n1, score_n2;       // Dimensions of scoring matrices
  long           score_size;               // Size of scoring matrices
  int            _flag;                    // Flag defining alignment mode
  AGEaligner     *_aux_aligner;            // Auxilary aligner

  // Coordiantes of the excised retion: left1, left2, right1, right2
  static const int MAX_BPOINTS = 100;
  int _bpoints[MAX_BPOINTS][4],_n_bpoints;
  int _max;                                // Maximum score

  short _match,_mismatch,_gap_open,_gap_extend; // Scoring parameters

public:
  AGEaligner(Sequence &s1,Sequence &s2);
  ~AGEaligner();
  bool align(Scorer &scr,int flag);
  void printAlignment();
  int score();

private:
  void printMatrix(unsigned short *matr);
  int findBPs(bool forward,int max_bps);
  int findInversionTDuplicationBPs(bool forward,int max_bps);
  int findIndelBPs(bool forward,int max_bps);
  Block *findLeftBlocks(int left1,int left2);
  Block *findRightBlocks(int right1,int right2);
  int  traceBackMaxima(int &lg1,int &lg2,int &rg1,int &rg2);
  int addBpoints(int lg1,int lg2,int rg1,int rg2); // 1 yes, 0 no, -1 overflow
  bool findAlignment(int bp_index = 0,bool forAlt = false);
  void calcScores(Scorer &scr);
  void calcMaxima();
  int calcWidth(int v1,int v2);
  int calcOutsideIdentity(string &seq,int bp1,int bp2);
  int calcInsideIdentity(string &seq,int bp1,int bp2);
  int calcIdentity(string &seq,int bp1,int bp2,int &left,int &right);
  void freeFragments(AliFragment *frags);
  bool getExcisedRange(AliFragment *f,int &s,int &e,
		       int add_s = 0,int add_e = 0);
};

class Block
{
private:
  Block *_next, *_prev;
  int    _start1,_start2,_length;
  string _desc;

public:
  Block(int st1,int st2,int l,string d = "") : _next(NULL), _prev(NULL),
					       _start1(st1), _start2(st2),
					       _length(l), _desc(d)
  {}

  string toString()
  {
    string ret = _desc;
    if (ret.length() > 0) ret += " ";
    ret += "block\t";
    stringstream tmp;
    tmp.str(""); tmp<<_start1;
    ret += tmp.str();  ret += "\t";
    tmp.str(""); tmp<<_start2;
    ret += tmp.str();  ret += "\t";
    tmp.str(""); tmp<<_length;
    ret += tmp.str();
    return ret;
  }

  Block *next(Block *b)
  {
    if (_next) return NULL;
    _next = b;
    b->prev(this);
    return _next;
  }

  Block *prev(Block *b)
  {
    if (_prev) return NULL;
    _prev = b;
    b->next(this);
    return _prev;
  }

  inline Block *next()   { return _next; }
  inline Block *prev()   { return _prev; }
  inline int    start1() { return _start1; }
  inline int    start2() { return _start2; }
  inline int    length() { return _length; }
};

class AliFragment
{
private:
  AliFragment *_next, *_prev;
  string       _ali1,  _ali2;
  int          _start1,_start2;
  int          _end1,  _end2;
public:
  AliFragment(string &ali1,string &ali2,
	      int start1,  int start2,
	      int end1,    int end2) : _next(NULL),     _prev(NULL),
				       _ali1(ali1),     _ali2(ali2),
				       _start1(start1), _start2(start2),
				       _end1(end1),     _end2(end2)
  {}

  AliFragment *next(AliFragment *b)
  {
    if (_next) return NULL;
    _next = b;
    b->prev(this);
    return _next;
  }

  AliFragment *prev(AliFragment *b)
  {
    if (_prev) return NULL;
    _prev = b;
    b->next(this);
    return _prev;
  }

  inline AliFragment *next()   { return _next; }
  inline AliFragment *prev()   { return _prev; }
  inline int          start1() { return _start1; }
  inline int          start2() { return _start2; }
  inline int          end1()   { return _end1; }
  inline int          end2()   { return _end2; }

  void countAligned(int &n_ali,int &n_id,int &n_gap)
  {
    n_ali = n_id = n_gap = 0;
    int len = _ali1.length();
    if (_ali2.length() < len) len = _ali2.length();
    for (int i = 0;i < len;i++) {
      n_ali++;
      if (Sequence::sameNuc(_ali1[i],_ali2[i]) > 0)               n_id++;
      if (Sequence::isGap(_ali1[i]) || Sequence::isGap(_ali2[i])) n_gap++;
    }
  }

  void printAlignment()
  {
    static const int WIDTH  = 60;
    static const int MARGIN =  9;

    int inc1 = 1; if (_end1 < _start1) inc1 = -1;
    int inc2 = 1; if (_end2 < _start2) inc2 = -1;

    string margin = "";
    for (int j = 0;j <= MARGIN;j++) margin += ' ';

    int n = _ali1.length(); if (_ali2.length() < n) n = _ali2.length();
    int ind1 = _start1, ind2 = _start2;
    for (int i = 0;i < n;i += WIDTH) {
      int st1 = ind1,st2 = ind2;
      string a1 = _ali1.substr(i,WIDTH);
      string a2 = _ali2.substr(i,WIDTH);
      int nuc1 = 0,nuc2 = 0;
      string match = "";
      for (int j = 0;j < a1.length();j++) {
	if (Sequence::sameNuc(a1[j],a2[j]) > 0) match += '|';
	else if (!Sequence::isGap(a1[j]) && 
		 !Sequence::isGap(a2[j]))       match += '.';
	else match += ' ';
	if (!Sequence::isGap(a1[j])) { nuc1++; ind1 += inc1; }
	if (!Sequence::isGap(a2[j])) { nuc2++; ind2 += inc2; }
      }

      cout<<endl;

      if (nuc1 > 0) {
	cout<<setw(MARGIN)<<st1<<' '<<a1;
	cout<<' '<<ind1 - inc1<<endl;
      } else cout<<margin<<a1<<endl;

      cout<<margin<<match<<endl;

      if (nuc2 > 0) {
	cout<<setw(MARGIN)<<st2<<' '<<a2;
	cout<<' '<<ind2 - inc2<<endl;
      } else cout<<margin<<a2<<endl;
    }
  }
};

#endif
