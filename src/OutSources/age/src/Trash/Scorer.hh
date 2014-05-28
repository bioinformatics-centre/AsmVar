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

#ifndef __SCORER_HH__
#define __SCORER_HH__

class Scorer
{
private:
  short _scores[256][256];
  short _match,_mismatch,_gap_open,_gap_extend;

public:
  Scorer(short match,short mismatch,short gap_open,short gap_extend) :
    _match(match),_mismatch(mismatch),
    _gap_open(gap_open),_gap_extend(gap_extend)
  {
    memset(_scores,0,256*256*sizeof(short));
    for (int i1 = 0;i1 < 256;i1++) {
      char c1 = (char)i1;
      if (Sequence::complement(c1) == c1) continue;
      for (int i2 = 0;i2 < 256;i2++) {
	char c2 = (char)i2;
	if (Sequence::complement(c2) == c2) continue;
	if (Sequence::sameNuc(c1,c2)) _scores[c1][c2] = match;
	else                          _scores[c1][c2] = mismatch;
      }
    }
  }
  
  inline short  getScore(char c1,char c2) { return _scores[c1][c2]; }
  inline short *getScores(char c)         { return &_scores[c][0]; }
  inline short  getMismatch()             { return _mismatch; }
  inline short  getMatch()                { return _match; }
  inline short  getGapOpen()              { return _gap_open; }
  inline short  getGapExtend()            { return _gap_extend; }
};

#endif
