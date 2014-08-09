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

#include "AGEaligner.h"

AGEaligner::AGEaligner(Sequence &s1,Sequence &s2) :
	_s1(s1), _s2(s2),
	_seq1(_s1.sequence()), _len1(_seq1.length()),
	_seq2(_s2.sequence()), _len2(_seq2.length()),
	_s1_rc(_s1.clone()),_seq1_rc(_s1_rc->sequence()),
	_frags(NULL), _frags_alt(NULL),_max(0), _flag(0), _aux_aligner(NULL),
	_match(0), _mismatch(0), _gap_open(0), _gap_extend(0),
	_n_bpoints(0), _is_set_align_result(false)
{
	// Constructing of reverse complement (RC) of first sequence 
	// RC could be needed for alignment in regions of inversion
	_s1_rc->revcom();

	score_n1 = _len1 + 2;
	score_n2 = _len2 + 2;
	score_size = (long)score_n1*(long)score_n2;

	// Scoring matrices
	f_score = new unsigned short[score_size];
	r_score = new unsigned short[score_size];
	memset(f_score,0,score_size*sizeof(unsigned short));
	memset(r_score,0,score_size*sizeof(unsigned short));

	// Trace matrix
	trace   = new unsigned char[score_size];
	memset(trace,0,score_size*sizeof(unsigned char));
}

AGEaligner::~AGEaligner()
{
	delete[] f_score;
	delete[] r_score;
	delete[] trace;
	delete   _s1_rc;
	freeFragments(_frags);
	freeFragments(_frags_alt);
	delete _aux_aligner;
}

void AGEaligner::freeFragments(AliFragment *frags)
{
	AliFragment *f = frags;
	while (f) {
		AliFragment *tmp = f;
		f = f->next();
		delete tmp;
	}
}

bool AGEaligner::align(Scorer &scr, int flag)
{
	_n_bpoints = 0; // Erase previous break points
	if (_len1 <= 0 || _len2 <= 0) return false;
	if (!(flag & ALL_FLAGS)) return false;

	// Deciding on whether need to use auxiliary aligner
	int aux_flag = 0;
	_flag = flag;
	if (flag & INVERSION_FLAG) {
		_flag    = INVL_FLAG;
		aux_flag = INVR_FLAG;
	}

#ifdef AGE_TIME
	timeval score_s,score_e;
	timeval max_s,  max_e;
	timeval find_s, find_e;
	timeval ali_s,  ali_e;

	gettimeofday(&score_s,NULL);
#endif

	_match      = scr.getMatch();
	_mismatch   = scr.getMismatch();
	_gap_open   = scr.getGapOpen();
	_gap_extend = scr.getGapExtend();

	calcScores(scr);
	//   printMatrix(f_score);
	//   printMatrix(r_score);

#ifdef AGE_TIME
	gettimeofday(&score_e,NULL);
	gettimeofday(&max_s,NULL);
#endif

	calcMaxima();
	//   printMatrix(f_score);
	//   printMatrix(r_score);

#ifdef AGE_TIME
	gettimeofday(&max_e,NULL);
	gettimeofday(&find_s,NULL);
#endif

	_max = findBPs(true,MAX_BPOINTS/2);
	_max = findBPs(false,MAX_BPOINTS /2);
	//cout<<"# bpoints = "<<_n_bpoints<<"\n";

#ifdef AGE_TIME
	gettimeofday(&find_e,NULL);
	gettimeofday(&ali_s,NULL);
#endif

	findAlignment(); // Finding alignment

#ifdef AGE_TIME
	gettimeofday(&ali_e,NULL);

	//   cout<<"Scoring  "<<score_e.tv_sec - score_s.tv_sec +
	//     (score_e.tv_usec - score_s.tv_usec)/1000000.<<" s"<<"\n";
	//   cout<<"Maxima   "<<max_e.tv_sec - max_s.tv_sec +
	//     (max_e.tv_usec - max_s.tv_usec)/1000000.<<" s"<<"\n";
	//   cout<<"Finding  "<<find_e.tv_sec - find_s.tv_sec +
	//     (find_e.tv_usec - find_s.tv_usec)/1000000.<<" s"<<"\n";
	//   cout<<"Aligning "<<ali_e.tv_sec - ali_s.tv_sec +
	//     (ali_e.tv_usec - ali_s.tv_usec)/1000000.<<" s"<<"\n";
#endif

	// Doing auxilary alignment 
	delete _aux_aligner;
	_aux_aligner = NULL;
	if (aux_flag) {
		_aux_aligner = new AGEaligner(_s1,_s2);
		if (!_aux_aligner->align(scr,aux_flag)) {
			delete _aux_aligner;
			_aux_aligner = NULL;
		}
	}

	return true;
}

int AGEaligner::score()
{
	if (_aux_aligner && _aux_aligner->score() > _max)
		return _aux_aligner->score();
	return _max;
}

bool AGEaligner::getExcisedRange(AliFragment *f,int &s,int &e,
		int add_s,int add_e)
{
	if (!f->next()) return false;
	int inc1 = 1; if (_s1.reverse()) inc1 = -1;
	if (_flag & INDEL_FLAG) {
		s = f->end1() + inc1;
		e = f->next()->start1() - inc1;
	} else if (_flag & TDUPLICATION_FLAG) {
		s = f->next()->start1();
		e = f->end1();
	} else if (_flag & INVL_FLAG) {
		s = f->end1() + inc1;
		e = f->next()->start1();
	} else if (_flag & INVR_FLAG) {
		s = f->end1();
		e = f->next()->start1() - inc1;
	} else return false;
	s += inc1*add_s;
	e += inc1*add_e;
	return true;
}

void AGEaligner::SetAlignResult() {

	_is_set_align_result = true;

	if (_aux_aligner && _aux_aligner->score() > _max) {
        _aux_aligner->printAlignment();
        return;
    }

	_align_result._score = score();

	int n_frg = 0;
    for (AliFragment *f = _frags; f; f = f->next()) n_frg++;
    int *n_ali = new int[n_frg + 1];
    int *n_id  = new int[n_frg + 1];
    int *n_gap = new int[n_frg + 1];
    n_ali[0] = n_id[0] = n_gap[0] = 0;
    int index = 1;
    for (AliFragment *f = _frags;f;f = f->next()) {
        f->countAligned(n_ali[index],n_id[index],n_gap[index]);
        n_ali[0] += n_ali[index];
        n_id[0]  += n_id[index];
        n_gap[0] += n_gap[index];
        index++;
    }

	_align_result._identity.clear();
    int identic = 0;
    if (n_ali[0] > 0) { 
		identic = (int)(100.*n_id[0]/n_ali[0] + 0.5);
		_align_result._identity.push_back(make_pair(n_id[0], identic));
	}
    if (n_frg > 1) {
        for (int i = 1;i < n_frg + 1; ++i) {
			_align_result._identity.push_back(
				make_pair(n_id[i], int(100.0*n_id[i]/n_ali[i] + 0.5)));
		}
    }

	_align_result._id1    = _s1.name(); // Id of first  Seq
	_align_result._id2    = _s2.name(); // Id of second Seq
	_align_result._strand = '.';
	_align_result._map.clear();
	_align_result._map_info.clear();
    // Printing aligned region coordinates
    if (n_ali[0] > 0) {

        for (AliFragment *f = _frags; f; f = f->next()) {

			// <First, second>
			MapData map1, map2;
			map1._id = _align_result._id1;
			map2._id = _align_result._id2;
			map1._start    = f->start1(); map1._end = f->end1();
			map2._start    = f->start2(); map2._end = f->end2();
			map1._sequence = f->ali1();
			map2._sequence = f->ali2();
			
			_align_result._map.push_back(make_pair(map1, map2));
			_align_result._map_info.push_back(f->mapInfo());
        }
		_align_result._strand = (_align_result._map[0].second._start > 
								 _align_result._map[0].second._end) ? '-' : '+';
    }

	int inc1 = 1; if (_s1.reverse()) inc1 = -1;
	int inc2 = 1; if (_s2.reverse()) inc2 = -1;
    if (n_frg > 1) {

		_align_result._is_alternative_align = (_n_bpoints > 1);

        // Printing sequence identity around breakpoints
		vector< pair<int,int> > ci_start1, ci_end1, ci_start2, ci_end2; 
		// Identity at breakpoints
        int left, right, s, e;
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {

            if (!getExcisedRange(f,s,e,0,1)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int homo_run = calcIdentity(_seq1,bp1,bp2,left,right);
			_align_result._homo_run_atbp1 = homo_run;

            if (homo_run > 0) {
				ci_start1.push_back(make_pair(s + inc1*left, s + inc1*(right - 1)));
				ci_end1.push_back  (make_pair(e + inc1*left, e + inc1*(right - 1)));
			}

            s = f->end2() + inc2, e = f->next()->start2();
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            homo_run = calcIdentity(_seq2,bp1,bp2,left,right);
			_align_result._homo_run_atbp2 = homo_run;

            if (homo_run > 0) {
				ci_start2.push_back(make_pair(s + inc2*left, s + inc2*(right - 1)));
				ci_end2.push_back  (make_pair(e + inc2*left, e + inc2*(right - 1)));
			}
        }

        // Identity outside breakpoints
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {

            if (!getExcisedRange(f,s,e,-1,1)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int homo_run = calcOutsideIdentity(_seq1,bp1,bp2);
			_align_result._homo_run_outbp1 = homo_run;

            if (homo_run > 0) {
				ci_start1.push_back(make_pair(s - inc1*(homo_run - 1), s));
				ci_end1.push_back  (make_pair(e, e + inc1*(homo_run - 1)));
			}

            s = f->end2(), e = f->next()->start2();
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            homo_run = calcOutsideIdentity(_seq2,bp1,bp2);
			_align_result._homo_run_outbp2 = homo_run; 
            if (homo_run > 0) {
				ci_start2.push_back(make_pair(s - inc2*(homo_run - 1), s));
				ci_end2.push_back  (make_pair(e, e + inc2*(homo_run - 1)));
			}
        }
		// Identity inside breakpoints
        for (AliFragment *f = _frags;f && f->next();f = f->next()) {
            if (!getExcisedRange(f,s,e)) continue;
            int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
            int homo_run = calcInsideIdentity(_seq1,bp1,bp2);
			_align_result._homo_run_inbp1 = homo_run;
            if (homo_run > 0) {
				ci_start1.push_back(make_pair(s, s + inc1*(homo_run - 1)));
				ci_end1.push_back  (make_pair(e - inc1*(homo_run - 1), e));
			}

            s = f->end2() + inc2, e = f->next()->start2() - inc2;
            bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
            homo_run = calcInsideIdentity(_seq2,bp1,bp2);
			_align_result._homo_run_inbp2 = homo_run;
            if (homo_run > 0) {
				ci_start2.push_back(make_pair(s, s + inc2*(homo_run - 1)));
				ci_end2.push_back  (make_pair(e - inc2*(homo_run - 1), e));
			}
        }

		_align_result._ci_start1 = Boundary(ci_start1);
		_align_result._ci_end1   = Boundary(ci_end1);
		_align_result._ci_start2 = Boundary(ci_start2);
		_align_result._ci_end2   = Boundary(ci_end2);
    }

    delete[] n_ali;
    delete[] n_id;
    delete[] n_gap;

	return;
}

void AGEaligner::printAlignment()
{
	if (_aux_aligner && _aux_aligner->score() > _max) {
		_aux_aligner->printAlignment();
		return;
	}

	const static string EXCISED_MESSAGE  = "EXCISED REGION";
	const static string EXCISED_MESSAGEs = "EXCISED REGION(S)";

	cout <<"\nMATCH = "   << _match    <<", ";
	cout <<"MISMATCH = "  << _mismatch <<", ";
	cout <<"GAP OPEN = "  << _gap_open <<", ";
	cout <<"GAP EXTEND = "<< _gap_extend;
	if      (_flag & INDEL_FLAG)        cout <<", INDEL";
	else if (_flag & INVL_FLAG)         cout <<", INVERSION";
	else if (_flag & INVR_FLAG)         cout <<", INVERSION";
	else if (_flag & TDUPLICATION_FLAG) cout <<", TDUPLICATION";
	cout << "\n\n";

	int inc1 = 1; if (_s1.reverse()) inc1 = -1;
	int inc2 = 1; if (_s2.reverse()) inc2 = -1;
	int s1 = _s1.start(),s2 =_s2.start();

	int e1 = s1 + inc1*(_seq1.length() - 1);
	int e2 = s2 + inc2*(_seq2.length() - 1);
	cout << "First  seq ["
		<< setw(calcWidth(s1,s2))  << s1 << ","
		<< setw(calcWidth(e1,e2))  << e1 << "] => "
		<< setw(9)<<_seq1.length() << " nucs '" << _s1.name() << "'\n";
	cout << "Second seq ["
		<< setw(calcWidth(s1,s2))  << s2 << ","
		<< setw(calcWidth(e1,e2))  << e2 << "] => "
		<< setw(9)<<_seq2.length() << " nucs '" << _s2.name() << "'\n\n";

	int n_frg = 0;
	for (AliFragment *f = _frags;f;f = f->next()) n_frg++;
	int *n_ali = new int[n_frg + 1];
	int *n_id  = new int[n_frg + 1];
	int *n_gap = new int[n_frg + 1];
	n_ali[0] = n_id[0] = n_gap[0] = 0;
	int index = 1;
	for (AliFragment *f = _frags;f;f = f->next()) {
		f->countAligned(n_ali[index],n_id[index],n_gap[index]);
		n_ali[0] += n_ali[index];
		n_id[0]  += n_id[index];
		n_gap[0] += n_gap[index];
		index++;
	}

	int identic = 0,gap = 0;
	if (n_ali[0] > 0) {
		identic = (int)(100.*n_id[0]/n_ali[0] + 0.5);
		gap     = (int)(100.*n_gap[0]/n_ali[0] + 0.5);
	}

	cout << "Score:   " << setw(9) << _max     << "\n";
	cout << "Aligned: " << setw(9) << n_ali[0] <<"        nucs\n";
	cout << "Identic: " << setw(9) << n_id[0]  <<" (" << setw(3) << identic << "%) nucs";
	if (n_frg > 1) {
		cout << " =>";
		for (int fr = 1;fr < n_frg + 1;fr++)
			cout << " " << setw(9) << n_id[fr] << " ("     << setw(3)
				<< (int)(100.*n_id[fr] / n_ali[fr] + 0.5) << "%)";
	}
	cout << "\nGaps:    "      << setw(9) << n_gap[0] << " (" << setw(3)
		<< gap << "%) nucs\n\n";

	// Printing aligned region coordinates
	if (n_ali[0] > 0) {
		cout << "Alignment:\n";
		cout << " first  seq =>  ";
		for (AliFragment *f = _frags;f;f = f->next()) {
			if (f->prev()) cout << " " << EXCISED_MESSAGE << " ";
			int ws = calcWidth(f->start1(),f->start2());
			int we = calcWidth(f->end1(),  f->end2());
			cout << "[" << setw(ws) << f->start1() << "," << setw(we) << f->end1() << "]";
		}
		cout << "\n second seq =>  ";
		for (AliFragment *f = _frags;f;f = f->next()) {
			if (f->prev()) cout << " " << EXCISED_MESSAGE << " ";
			int ws = calcWidth(f->start1(),f->start2());
			int we = calcWidth(f->end1(),  f->end2());
			cout << "[" << setw(ws) << f->start2() << "," << setw(we) << f->end2() << "]";
		}
		cout << "\n";
	}

	int s,e,len;
	if (n_frg > 1) {
		cout << "\n" << EXCISED_MESSAGEs << ":\n";
		for (AliFragment *f = _frags;f && f->next();f = f->next()) {
			if (!getExcisedRange(f,s,e)) continue;
			cout << " first  seq => ";
			len = abs(s - e - inc1);
			cout << setw(9) << len << " nucs";
			if (len > 0) cout << " [" << s << "," << e << "]";
			cout << "\n";
			cout << " second seq => ";
			s = f->end2() + inc2;
			e = f->next()->start2() - inc2;
			len = abs(s - e - inc2);
			cout << setw(9) << len << " nucs";
			if (len > 0) cout << " [" << s << "," << e << "]";
			cout << "\n";
		}

		if (_n_bpoints > 1) cout<<"ALTERNATIVE REGION(S): "<<_n_bpoints - 1<<"\n";
		for (int i = _n_bpoints - 1;i >= 1;i--) {
			if (!findAlignment(i,true)) continue;  // Fiding alternative alignment
			if (!_frags_alt) continue;
			for (AliFragment *f = _frags_alt;f && f->next();f = f->next()) {
				if (!getExcisedRange(f,s,e)) continue;
				cout<<" first  seq => ";
				len = abs(s - e - inc1);
				cout<<setw(9)<<len<<" nucs";
				if (len > 0) cout<<" ["<<s<<","<<e<<"]";
				cout<<"\n";
				cout<<" second seq => ";
				s = f->end2() + inc2;
				e = f->next()->start2() - inc2;
				len = abs(s - e - inc2);
				cout<<setw(9)<<len<<" nucs";
				if (len > 0) cout<<" ["<<s<<","<<e<<"]";
				cout<<"\n";
			}
			break;
		}

		// Printing sequence identity around breakpoints
		cout<<"\n"<<"Identity at breakpoints: "<<"\n";
		int left,right,s,e;
		for (AliFragment *f = _frags;f && f->next();f = f->next()) {
			if (!getExcisedRange(f,s,e,0,1)) continue;
			int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
			int homo_run = calcIdentity(_seq1,bp1,bp2,left,right);
			cout<<" first  seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s + inc1 * left<<","<<s + inc1*(right - 1)<<"] to"
					<<" ["<<e + inc1*left<<","<<e + inc1*(right - 1)<<"]";
			cout<<"\n";
			s = f->end2() + inc2, e = f->next()->start2();
			bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
			homo_run = calcIdentity(_seq2,bp1,bp2,left,right);
			cout<<" second seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s + inc2*left<<","<<s + inc2*(right - 1)<<"] to"
					<<" ["<<e + inc2*left<<","<<e + inc2*(right - 1)<<"]";
			cout<<"\n";
		}

		cout<<"Identity outside breakpoints: "<<"\n";
		for (AliFragment *f = _frags;f && f->next();f = f->next()) {
			if (!getExcisedRange(f,s,e,-1,1)) continue;
			int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
			int homo_run = calcOutsideIdentity(_seq1,bp1,bp2);
			cout<<" first  seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s - inc1*(homo_run - 1)<<","<<s<<"] to"
					<<" ["<<e<<","<<e + inc1*(homo_run - 1)<<"]";
			cout<<"\n";
			s = f->end2(), e = f->next()->start2();
			bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
			homo_run = calcOutsideIdentity(_seq2,bp1,bp2);
			cout<<" second seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s - inc2*(homo_run - 1)<<","<<s<<"] to"
					<<" ["<<e<<","<<e + inc2*(homo_run - 1)<<"]";
			cout<<"\n";
		}

		cout<<"Identity inside breakpoints: "<<"\n";
		for (AliFragment *f = _frags;f && f->next();f = f->next()) {
			if (!getExcisedRange(f,s,e)) continue;
			int bp1 = inc1*(s - _s1.start()), bp2 = inc1*(e - _s1.start());
			int homo_run = calcInsideIdentity(_seq1,bp1,bp2);
			cout<<" first  seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s<<","<<s + inc1*(homo_run - 1)<<"] to"
					<<" ["<<e - inc1*(homo_run - 1)<<","<<e<<"]";
			cout<<"\n";
			s = f->end2() + inc2, e = f->next()->start2() - inc2;
			bp1 = inc2*(s - _s2.start()), bp2 = inc2*(e - _s2.start());
			homo_run = calcInsideIdentity(_seq2,bp1,bp2);
			cout<<" second seq => "<<setw(9)<<homo_run<<" nucs";
			if (homo_run > 0)
				cout<<" ["<<s<<","<<s + inc2*(homo_run - 1)<<"] to"
					<<" ["<<e - inc2*(homo_run - 1)<<","<<e<<"]";
			cout<<"\n";
		}
	}

	// Printing actual alignment
	for (AliFragment *f = _frags;f;f = f->next()) {
		if (f != _frags) cout<<"\n"<<EXCISED_MESSAGE<<"\n";
		f->printAlignment();
	}

	delete[] n_ali;
	delete[] n_id;
	delete[] n_gap;
}

int AGEaligner::calcWidth(int v1,int v2)
{
	int width1 = 0; if (v1 < 0) v1++;
	while (v1 != 0) { width1++; v1 /= 10; }
	int width2 = 0; if (v2 < 0) v2++;
	while (v2 != 0) { width2++; v2 /= 10; }
	if (width1 > width2) return width1;
	return width2;
}

int AGEaligner::calcOutsideIdentity(string &seq,int bp1,int bp2)
{
	// bp1 and bp2 are zero based
	if (bp1 >= bp2) return 0;
	int n = seq.length(),d1 = bp1 + 1, d2 = n - bp2;
	int n_check = d1; if (d2 < d1) n_check = d2;
	int start1 = bp1 - n_check + 1;
	for (int i = 0;i < n_check;i++) {
		int n_same = 0,start2 = bp2 - i;
		for (int j = i;j < n_check;j++)
			if (Sequence::sameNuc(seq[start1 + j],seq[start2 + j])) n_same++;
			else break;
		if (n_same == n_check - i) return n_same;
	}

	return 0;
}

int AGEaligner::calcInsideIdentity(string &seq, int bp1, int bp2)
{
	// bp1 and bp2 are zero based
	if (bp1 >= bp2) return 0;
	int n_check = (bp2 - bp1 + 1)/2, start2 = bp2 - n_check + 1;
	for (int i = 0;i < n_check;i++) {
		int n_same = 0,start1 = bp1 - i;
		for (int j = i;j < n_check;j++) {
			if (Sequence::sameNuc(seq[start1 + j], seq[start2 + j])) n_same++;
			else break;
		}
		if (n_same == n_check - i) return n_same;
	}

	return 0;
}

int AGEaligner::calcIdentity(string &seq,int bp1,int bp2,int &left,int &right)
{
	// bp1 and bp2 are zero based
	if (bp1 >= bp2) return 0;
	int n = seq.length(),delta = bp2 - bp1;
	if (delta == 0) return 0;
	int start = bp1 - 1;
	left = right = 0;
	while (start >= 0 && Sequence::sameNuc(seq[start],seq[start + delta])) {
		start--;
		left--;
	}

	int end = bp1;
	while (end < n && Sequence::sameNuc(seq[end],seq[end + delta])) {
		end++;
		right++;
	}

	return right - left;
}

pair<int, int> AGEaligner::Boundary(vector< pair<int, int> > reg) {

    pair<int, int> boundary(make_pair(0, 0));
	if (reg.size() > 0) {

		boundary = reg[0];
    	for (size_t i(1); i < reg.size(); ++i) {

			if (boundary.first  >  reg[i].first) boundary.first  = reg[i].first;
			if (boundary.second < reg[i].second) boundary.second = reg[i].second;
    	}
	}

    return boundary;
}

void AGEaligner::printMatrix(unsigned short *matr)
{
	long pos = 0;
	for (int i2 = 0;i2 < score_n2;i2++) {
		for (int i1 = 0;i1 < score_n1;i1++)
			cout<<setw(4)<<matr[pos++]<<" ";
		cout<<"\n";
	}
	cout<<"\n";
}

int AGEaligner::findBPs(bool forward,int max_bps)
{
	int max = 0;
	if      (_flag & INDEL_FLAG)
		max = findIndelBPs(forward,max_bps);
	else if (_flag & (INVL_FLAG | INVR_FLAG | TDUPLICATION_FLAG))
		max = findInversionTDuplicationBPs(forward,max_bps);
	else return max;

	return max;
}

int AGEaligner::findInversionTDuplicationBPs(bool forward,int max_bps)
{
	int max = 0,n_add = 0;
	unsigned short *f_maxima = f_score, *r_maxima = r_score;
	if (forward) {
		long posf = score_n1 - 2,posr = score_n1 + 1;
		for (int i2 = 0;i2 < score_n2 - 1;i2++) {
			unsigned short val = f_maxima[posf] + r_maxima[posr];
			if (val > max) max = val;
			posf += score_n1;
			posr += score_n1;
		}
		posf = score_n1 - 2,posr = score_n1 + 1;
		for (int i2 = 0;i2 < score_n2 - 1;i2++) {
			if (f_maxima[posf] + r_maxima[posr] == max) 
				for (int j1 = 0;j1 < score_n1 - 1;j1++) {
					int ind1 = posr - score_n1 + j1 - 1;
					if ((trace[ind1] & FM_MASK) != 0) continue;
					for (int j2 = 0;j2 < score_n1 - 1;j2++) {
						int ind2 = posr + j2;
						if (f_maxima[ind1] + r_maxima[ind2] == max) {
							int lg1 = j1,    lg2 = i2;
							int rg1 = j2 + 1,rg2 = i2 + 1;
							traceBackMaxima(lg1,lg2,rg1,rg2);
							int res = addBpoints(lg1,lg2,rg1,rg2);
							if (res > 0) n_add++;
							if (n_add >= max_bps || res < 0) return max;
						}
					}
				}
			posf += score_n1;
			posr += score_n1;
		}
	} else { // Reverse
		long posf = score_size - score_n1 - 2,posr = score_size - score_n1 + 1;
		for (int i2 = score_n2 - 1;i2 > 0;i2--) {
			unsigned short val = f_maxima[posf] + r_maxima[posr];
			if (val > max) max = val;
			posf -= score_n1;
			posr -= score_n1;
		}

		posf = score_size - score_n1 - 2,posr = score_size - score_n1 + 1;
		for (int i2 = score_n2 - 1;i2 > 0;i2--) {
			if (f_maxima[posf] + r_maxima[posr] == max)
				for (int j2 = score_n1 - 1;j2 > 0;j2--) {
					int ind2 = posr + j2 - 1;
					if ((trace[ind2] & RM_MASK) != 0) continue;
					for (int j1 = score_n1 - 1;j1 > 0;j1--) {
						int ind1 = posr - score_n1 + j1 - 2;
						if (f_maxima[ind1] + r_maxima[ind2] == max) {
							int lg1 = j1 - 1,lg2 = i2 - 1;
							int rg1 = j2,    rg2 = i2;
							traceBackMaxima(lg1,lg2,rg1,rg2);
							int res = addBpoints(lg1,lg2,rg1,rg2);
							if (res > 0) n_add++;
							if (n_add >= max_bps || res < 0) return max;
						}
					}
				}
			posf -= score_n1;
			posr -= score_n1;
		}
	}

	return max;
}

int AGEaligner::findIndelBPs(bool forward,int max_bps)
{
	int max = 0,n_add = 0;
	unsigned short *f_maxima = f_score, *r_maxima = r_score;
	if (forward) {
		unsigned short *maxf = f_score,*maxr = &r_score[score_n1 + 1];
		for (int i2 = 0;i2 < score_n2 - 1;i2++) {
			for (int i1 = 0;i1 < score_n1 - 1;i1++) {
				unsigned short val = *maxf + *maxr;
				if (val > max) max = val;
				maxf++;
				maxr++;
			}
			maxf++;
			maxr++;
		}
		maxf = f_score, maxr = &r_score[score_n1 + 1];
		unsigned char *tr = trace;
		for (int i2 = 0;i2 < score_n2 - 1;i2++) {
			for (int i1 = 0;i1 < score_n1 - 1;i1++) {
				if ((*tr & FM_MASK) == 0 && *maxf + *maxr == max) {
					int lg1 = i1,     lg2 = i2;
					int rg1 = i1 + 1, rg2 = i2 + 1;
					traceBackMaxima(lg1,lg2,rg1,rg2);
					int res = addBpoints(lg1,lg2,rg1,rg2);
					if (res > 0) n_add++;
					if (n_add >= max_bps || res < 0) return max;
				}
				maxf++;
				maxr++;
				tr++;
			}
			maxf++;
			maxr++;
			tr++;
		}
	} else { // Reverse
		unsigned short *maxf = &f_score[score_size - score_n1 - 2];
		unsigned short *maxr = &r_score[score_size - 1];
		for (int i2 = score_n2 - 1;i2 > 0;i2--) {
			for (int i1 = score_n1 - 1;i1 > 0;i1--) {
				unsigned short val = *maxf + *maxr;
				if (val > max) max = val;
				maxf--;
				maxr--;
			}
			maxf--;
			maxr--;
		}
		maxf = &f_score[score_size - score_n1 - 2];
		maxr = &r_score[score_size - 1];
		unsigned char *tr = &trace[score_size - 1];
		for (int i2 = score_n2 - 1;i2 > 0;i2--) {
			for (int i1 = score_n1 - 1;i1 > 0;i1--) {
				if ((*tr & RM_MASK) == 0 && *maxf + *maxr == max) {
					int lg1 = i1 - 1, lg2 = i2 - 1;
					int rg1 = i1,     rg2 = i2;
					traceBackMaxima(lg1,lg2,rg1,rg2);
					int res = addBpoints(lg1,lg2,rg1,rg2);
					if (res > 0) n_add++;
					if (n_add >= max_bps || res < 0) return max;
				}
				maxf--;
				maxr--;
				tr--;
			}
			maxf--;
			maxr--;
			tr--;
		}
	}

	return max;
}

int AGEaligner::traceBackMaxima(int &lg1,int &lg2,int &rg1,int &rg2)
{
	long pos = (long)lg2*(long)score_n1 + (long)lg1;
	unsigned char tr = trace[pos] & FM_MASK;
	while (tr != 0) {
		if (tr == FM_VERTICAL) {
			pos -= score_n1;
			lg2--;
		} else if (tr == FM_HORIZONTAL) {
			pos--;
			lg1--;
		} else cerr<<"Internal error (1) in tracking."<<"\n";
		tr = trace[pos] & FM_MASK;
	}

	pos = (long)rg2*(long)score_n1 + (long)rg1;
	tr = trace[pos] & RM_MASK;
	while (tr != 0) {
		if (tr == RM_VERTICAL) {
			pos += score_n1;
			rg2++;
		} else if (tr == RM_HORIZONTAL) {
			pos++;
			rg1++;
		} else cerr<<"Internal error (2) in tracking."<<"\n";
		tr = trace[pos] & RM_MASK;
	}
}

int AGEaligner::addBpoints(int lg1,int lg2,int rg1,int rg2)
{
	//cout<<lg1<<" "<<lg2<<" "<<rg1<<" "<<rg2<<"\n";
	if (lg1 == 0 || lg2 == 0) lg1 = lg2 = 0;
	if (rg1 == score_n1 - 1 || lg2 == score_n2 - 1) {
		rg1 = score_n1 - 1;
		rg2 = score_n2 - 1;
	}
	for (int i = 0;i < _n_bpoints;i++)
		if (_bpoints[i][0] - lg1 == _bpoints[i][2] - rg1 &&
				_bpoints[i][1] - lg2 == _bpoints[i][3] - rg2)
			return 0;

	if (_n_bpoints < MAX_BPOINTS) {
		_bpoints[_n_bpoints][0] = lg1;
		_bpoints[_n_bpoints][1] = lg2;
		_bpoints[_n_bpoints][2] = rg1;
		_bpoints[_n_bpoints][3] = rg2;
		_n_bpoints++;
		return 1;
	}

	//cerr<<"Buffer for breakpoints exceeded. Skipping ..."<<"\n";
	return -1;
}

bool AGEaligner::findAlignment(int bp_index,bool forAlt)
{
	if (forAlt) {
		freeFragments(_frags_alt);
		_frags_alt = NULL;
	} else {
		freeFragments(_frags);
		_frags = NULL;
	}

	if (bp_index > _n_bpoints) return false;

	static const int MAX_N_FRAGS = 2;
	int lg1 = _bpoints[bp_index][0];
	int lg2 = _bpoints[bp_index][1];
	int rg1 = _bpoints[bp_index][2];
	int rg2 = _bpoints[bp_index][3];
	Block *bs[MAX_N_FRAGS] = {findLeftBlocks(lg1,lg2),findRightBlocks(rg1,rg2)};

	int start1 = 0,start2 = 0;
	int end1   = 0,end2   = 0;

	AliFragment *frags = NULL;
	string ali1(""),ali2("");
	for (int b = 0;b < MAX_N_FRAGS;b++) {
		Block *bls = bs[b]; if (!bls) continue;
		bool use_rc1 = ((_flag & INVL_FLAG) && b == 1) ||
			((_flag & INVR_FLAG) && b == 0);
		Sequence *s1 = &_s1,  *s2   = &_s2;
		string *seq1 = &_seq1,*seq2 = &_seq2;
		if (use_rc1) { seq1 = &_seq1_rc; s1 = _s1_rc; }
		int inc1 = 1; if (s1->reverse()) inc1 = -1;
		int inc2 = 1; if (s2->reverse()) inc2 = -1;
		start1 = s1->start() + inc1*(bls->start1() - 1);
		start2 = s2->start() + inc2*(bls->start2() - 1);
		int nuc_ind = 2*b;
		int ind1 = bls->start1() - 1,ind2 = bls->start2() - 1;
		while (bls) {
			end1 = s1->start() + inc1*(bls->start1() + bls->length() - 2);
			end2 = s2->start() + inc2*(bls->start2() + bls->length() - 2);
			int st1 = bls->start1() - 1, st2 = bls->start2() - 1;
			int len = bls->length();
			while (ind1 < st1) {
				ali1 += (*seq1)[ind1++];
				ali2 += Sequence::gap();
			}
			while (ind2 < st2) {
				ali1 += Sequence::gap();
				ali2 += (*seq2)[ind2++];
			}
			for (int i = 0;i < len;i++) {
				ali1 += (*seq1)[ind1++];
				ali2 += (*seq2)[ind2++];
			}
			bls = bls->next();
		}
		if (ali1.length() > 0 && ali2.length() > 0 &&
				ali1.length() == ali2.length()) {
			AliFragment *f = new AliFragment(ali1,ali2,start1,start2,end1,end2);
			if (!frags) frags = f;
			else         frags = frags->next(f);
		}
		ali1 = ali2 = "";
	}

	while (frags && frags->prev()) frags = frags->prev();

	// Deleting blocks
	for (int i = 0;i < MAX_N_FRAGS;i++) {
		Block *bls = bs[i];
		while (bls) {
			Block *tmp = bls;
			bls = bls->next();
			delete tmp;
		}
	}
	if (forAlt) _frags_alt = frags;
	else        _frags     = frags;

	return true;
}

Block *AGEaligner::findLeftBlocks(int left1,int left2)
{
	Block *bls = NULL;

	// Producing alignment blocks for left alignment
	int n_save = _len1; if (_len2 > n_save) n_save = _len2;
	int *p1 = new int[n_save];
	int *p2 = new int[n_save];
	n_save = 0;
	long pos = (long)left2*(long)score_n1 + (long)left1;
	unsigned char tr = trace[pos] & FS_MASK;
	while (tr != 0) {
		if (tr == FS_DIAGONAL) {
			p1[n_save] = left1--;
			p2[n_save] = left2--;
			n_save++;
			pos -= (score_n1 + 1);
		} else if (tr == FS_HORIZONTAL) {
			left1--;
			pos--;
		} else if (tr == FS_VERTICAL) {
			left2--;
			pos -= score_n1;
		} else cerr<<"Internal error (3) in tracking."<<"\n";
		tr = trace[pos] & FS_MASK;
	}

	int start1,start2;
	if (n_save > 0) {
		int index = n_save - 1,n = 1;
		start1 = p1[index];
		start2 = p2[index];
		index--;
		while (index >= 0) {
			if (p1[index] - p1[index + 1] == 1 &&
					p2[index] - p2[index + 1] == 1) n++;
			else {
				Block *b = new Block(start1,start2,n,"left");
				if (bls == NULL) bls = b;
				else             bls = bls->next(b);
				start1 = p1[index];
				start2 = p2[index];
				n = 1;
			}
			index--;
		}
		Block *b = new Block(start1,start2,n,"left");
		if (bls == NULL) bls = b;
		else             bls = bls->next(b);
	}

	delete[] p1;
	delete[] p2;

	if (bls) while (bls->prev()) bls = bls->prev(); // Looping to the first
	return bls;
}

Block *AGEaligner::findRightBlocks(int right1,int right2)
{
	Block *bls = NULL;

	// Producing alignment blocks for right alignment
	long pos = (long)right2*(long)score_n1 + (long)right1;
	int start1,start2,n = 0;
	unsigned char tr = trace[pos] & RS_MASK;
	while (tr != 0) {
		if (tr == RS_DIAGONAL) {
			if (right1 <= _len1 && right2 <= _len2) {
				if (n == 0) {
					start1 = right1;
					start2 = right2;
					n = 1;
				} else n++;
			}
			right1++;
			right2++;
			pos += (score_n1 + 1);
		} else if (tr == RS_HORIZONTAL) {
			right1++;
			pos++;
			if (n > 0) {
				Block *b = new Block(start1,start2,n,"right");
				if (bls == NULL) bls = b;
				else             bls = bls->next(b);
				n = 0;
			}
		} else if (tr == RS_VERTICAL) {
			right2++;
			pos += score_n1;
			if (n > 0) {
				Block *b = new Block(start1,start2,n,"right");
				if (bls == NULL) bls = b;
				else             bls = bls->next(b);
				n = 0;
			}
		} else cerr<<"Internal error (4) in tracking."<<"\n";
		tr = trace[pos] & RS_MASK;
	}

	if (n > 0) {
		Block *b = new Block(start1,start2,n,"right");
		if (bls == NULL) bls = b;
		else             bls = bls->next(b);
	}

	if (bls) while (bls->prev()) bls = bls->prev();
	return bls;
}

void AGEaligner::calcMaxima()
{
	int thread1_i2 = -100,thread2_i2 = -100;
#ifdef OMP
#pragma omp shared(thread1_i2,thread2_i2)
	omp_set_num_threads(2);
#pragma omp parallel sections
	{
#pragma omp section
#endif

		{
			unsigned char  *tr = &trace[1];
			for (int i = 1;i < score_n1;i++) {
				*tr |= FM_HORIZONTAL;
				tr++;
			}
			for (int i = 1;i < score_n2 - 1;i++) {
				*tr |= FM_VERTICAL;
				tr += score_n1;
			}

			// Assigning maxima for leading matrices (forward direction)
			tr = &trace[score_n1 + 1];
			unsigned short *f  = &f_score[score_n1 + 1];
			unsigned short *f1 = &f_score[score_n1];
			unsigned short *f2 = &f_score[1];
			for (int i2 = 0;i2 < _len2;i2++) {
				thread1_i2 = i2;           // Indicate the row is in use
#ifdef OMP
#pragma omp flush(thread1_i2,thread2_i2)
#endif
				while (i2 == thread2_i2) ; // Wait for other thread
				for (int i1 = 0;i1 < _len1;i1++) {
					unsigned char  curr_trace = 0;
					if (*f2 > *f && i2 > 0) {
						*f = *f2;
						curr_trace = FM_VERTICAL;
					}
					if (*f1 > *f && i1 > 0) {
						*f = *f1;
						curr_trace = FM_HORIZONTAL;
					}
					*tr |= curr_trace;
					f++;
					f2++;
					f1++;
					tr++;
				}
				f++;  f++;
				f2++; f2++;
				f1++; f1++;
				tr++; tr++;
			}
			thread1_i2 = -100;
		}

#ifdef OMP
#pragma omp section
#endif

		{
			unsigned char  *tr = &trace[score_size - 2];
			for (int i = 1;i < score_n1;i++) {
				*tr |= FM_HORIZONTAL;
				tr--;
			}
			for (int i = 1;i < score_n2 - 1;i++) {
				*tr |= FM_VERTICAL;
				tr -= score_n1;
			}

			// Assigning maxima for trailing matrices (reverse direction)
			tr = &trace[score_size - score_n1 - 2];
			unsigned short *f  = &r_score[score_size - score_n1 - 2];
			unsigned short *f1 = &r_score[score_size - score_n1 - 1];
			unsigned short *f2 = &r_score[score_size - 2];
			int len1_tmp = _len1 - 1,len2_tmp = _len2 - 1;
			for (int i2 = len2_tmp;i2 >= 0;i2--) {
				thread2_i2 = i2;           // Indicate the row is in use
#ifdef OMP
#pragma omp flush(thread1_i2,thread2_i2)
#endif
				long tmp = 0;
				while (i2 == thread1_i2 &&
						tmp < score_size) tmp++; // Wait for other thread
				for (int i1 = len1_tmp;i1 >= 0;i1--) {
					unsigned char  curr_trace = 0;
					if (*f2 > *f && i2 < len2_tmp) {
						*f = *f2;
						curr_trace = RM_VERTICAL;
					}
					if (*f1 > *f && i1 < len1_tmp) {
						*f = *f1;
						curr_trace = RM_HORIZONTAL;
					}
					*tr |= curr_trace;
					f--;
					f2--;
					f1--;
					tr--;
				}
				f--;  f--;
				f2--; f2--;
				f1--; f1--;
				tr--; tr--;
			}
			thread2_i2 = -100;
		}

#ifdef OMP
	}
#endif

} // calcMaxima(...)

void AGEaligner::calcScores(Scorer &scr)
{
	int thread1_i2 = -100,thread2_i2 = -100;
#ifdef OMP
#pragma omp shared(thread1_i2,thread2_i2);
	omp_set_num_threads(2);
#pragma omp parallel sections
	{
#pragma omp section
#endif
		{
			// Filling in the forward direction
			unsigned short *fs1 = &f_score[score_n1], *fs2 = f_score;
			unsigned char  *tr1 = &trace[score_n1],   *tr2 = &trace[1];
			const char *c2 = _seq2.c_str();
			for (int i2 = 0;i2 < _len2;i2++) {
				thread1_i2 = i2;           // Indicate the row is in use
#ifdef OMP
#pragma omp flush(thread1_i2,thread2_i2)
#endif
				while (i2 == thread2_i2) ; // Wait for other thread
				const char *c1 = _seq1.c_str();
				if (_flag & INVR_FLAG) c1 = _seq1_rc.c_str();
				short *scs = scr.getScores(*c2);
				for (int i1 = 0;i1 < _len1;i1++) {
					int score = *(fs2++) + scs[*c1++];
					int sg2 = *fs2;
					unsigned char tr = *(tr2++) & FS_MASK;
					if (tr == FS_HORIZONTAL ||tr == FS_VERTICAL)
						sg2    += scr.getGapExtend();
					else sg2 += scr.getGapOpen();
					int sg1 = *fs1++;
					tr = *(tr1++) & FS_MASK;
					if (tr == FS_HORIZONTAL || tr == FS_VERTICAL)
						sg1    += scr.getGapExtend();
					else sg1 += scr.getGapOpen();
					unsigned char curr_trace = FS_DIAGONAL;
					if (sg1 >= score) { score = sg1; curr_trace = FS_HORIZONTAL; }
					if (sg2 >= score) { score = sg2; curr_trace = FS_VERTICAL; }
					if (score < 0)    { score = 0;   curr_trace = 0; }
					*tr1 |= curr_trace;
					*fs1  = (unsigned short)score;
				}
				fs1++; fs1++;
				fs2++; fs2++;
				tr1++; tr1++;
				tr2++; tr2++;
				c2++;
			}
			thread1_i2 = -100;
		}

#ifdef OMP
#pragma omp section
#endif

		{
			// Filling in the reverse direction
			unsigned short *rs1 = &r_score[score_size - 1 - score_n1];
			unsigned short *rs2 = &r_score[score_size - 1];
			unsigned char  *tr1 = &trace[score_size - 1 - score_n1];
			unsigned char  *tr2 = &trace[score_size - 2];
			const char *c2 = &(_seq2.c_str()[_len2 - 1]);
			for (int i2 = _len2 - 1;i2 >= 0;i2--) {
				thread2_i2 = i2;           // Indicate the row is in use
#ifdef OMP
#pragma omp flush(thread1_i2,thread2_i2)
#endif
				long tmp = 0;
				while (i2 == thread1_i2 &&
						tmp < score_size) tmp++; // Wait for other thread
				const char *c1 = &(_seq1.c_str()[_len1 - 1]);
				if (_flag & INVL_FLAG) c1 = &(_seq1_rc.c_str())[_len1 - 1];
				short *scs = scr.getScores(*c2);
				for (int i1 = _len1 - 1;i1 >= 0;i1--) {
					int score = *(rs2--) + scs[*c1--];
					int sg2 = *rs2;
					unsigned char tr = *(tr2--) & RS_MASK;
					if (tr == RS_HORIZONTAL || tr == RS_VERTICAL)
						sg2    += scr.getGapExtend();
					else sg2 += scr.getGapOpen();
					int sg1 = *rs1--;
					tr = *(tr1--) & RS_MASK;
					if (tr == RS_HORIZONTAL || tr == RS_VERTICAL)
						sg1    += scr.getGapExtend();
					else sg1 += scr.getGapOpen();
					unsigned char curr_trace = RS_DIAGONAL;
					if (sg1 >= score) { score = sg1; curr_trace = RS_HORIZONTAL; }
					if (sg2 >= score) { score = sg2; curr_trace = RS_VERTICAL; }
					if (score < 0)    { score = 0;   curr_trace = 0; }
					*tr1 |= curr_trace;
					*rs1 =  (unsigned short)score;
				}
				rs1--; rs1--;
				rs2--; rs2--;
				tr1--; tr1--;
				tr2--; tr2--;
				c2--;
			}
			thread2_i2 = -100;
		}

#ifdef OMP
	}
#endif
} // calcScores(...)

