/**********************************
 * Author : Shujia Huang
 * Date   : 2014-08-05 22:49:56
 **********************************/

#include "AgeAlignment.h"

AgeAlignment::AgeAlignment() : isInit(false) {}

bool AgeAlignment::Init(VarUnit &v, AgeOption op){
	vu     = v;
	opt    = op;
	isInit = true;
}

void ExtendVariant(unsigned long int tarFaSize, 
				   unsigned long int qryFaSize,
				   int extandFlankSzie){
	
	if (!isInit) { 
		cerr << "[ERROR] You should init AgeAlignment before calling Align()\n";
		exit(1);
	}
	
	vu.target.start -= extandFlankSzie;
	vu.target.end   += extandFlankSzie;
	vu.query.start  -= extandFlankSzie;
	vu.query.end    += extandFlankSzie;

	if (vu.target.start < 1) vu.target.start = 1;
    if (vu.target.end   > tarFaSize) vu.target.end = tarFaSize;	
    if (vu.query.start  < 1) vu.query.start = 1; 
    if (vu.query.end    > qryFaSize) vu.query.end  = qryFaSize;

	return;
}

bool AgeAlignment::Align(string &tarFa, string &qryFa){

	if (!isInit) { 
		cerr << "[ERROR] You should init AgeAlignment before calling Align()\n";
		exit(1);
	}
	ExtendVariant(tarFa.length(), qryFa.length(), opt.extendVarFlankSzie);

    // Do not AGE if the memory cost is bigger than 10G.
    if (IsHugeMemory(vu.target.end - vu.target.start, vu.query.end - vu.query.start))
		return false;

	Sequence *tarSeq = new Sequence(tarFa, vu.target.id, 1, false);
	Sequence *qrySeq = new Sequence(qryFa, vu.query.id , 1, false);
	Sequence *tar = tarSeq->substr(vu.target.start, vu.target.end);
	Sequence *qry = qrySeq->substr(vu.query.start , vu.query.end );

#ifdef AGE_TIME
timeval ali_s,ali_e;
gettimeofday(&ali_s, NULL);
#endif

	bool isalign(true);
	Scorer scr(opt.match, opt.mismatch, opt.gapOpen, opt.gapExtend);
	if (opt.both) {
            
        Sequence *qryClone = qry->clone(); 
        qryClone->revcom();
        AGEaligner aligner1(*tar, *qry);
        AGEaligner aligner2(*tar, *qryClone);                                         
        bool res1 = aligner1.align(scr, flag);
		bool res2 = aligner2.align(scr, flag);

		if (!res1 && !res2) {
			cerr<<"No alignment made.\n";
			isalign = false;
		} else if (aligner1.score() >= aligner2.score()) {
			aligner1.printAlignment();
        } else {
			aligner2.printAlignment();
        } 
        delete qryClone;

    } else {

        AGEaligner aligner(*tar, *qry);
        if (aligner.align(scr, flag)){
            aligner.printAlignment();
        } else {
            cerr << "No alignment made.\n";
			isalign = false; 
        }
    }

#ifdef AGE_TIME
gettimeofday(&ali_e, NULL);
cout << "\nAlignment time is " << ali_e.tv_sec - ali_s.tv_sec
    + (ali_e.tv_usec - ali_s.tv_usec)/1e+6 <<" s\n\n";
#endif

	Sequence::deleteSequences(tarSeq);
	Sequence::deleteSequences(qrySeq);

	return isalign;

}

bool IsHugeMemory(unsigned long int n, unsigned long int m){ 

	return (5 * n * m / 1000000000 > 10); // 10G
}

