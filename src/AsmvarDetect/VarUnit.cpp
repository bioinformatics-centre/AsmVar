#include "VarUnit.h"

VarUnit::VarUnit(){

	target.id = "-"; 
	query.id  = "-"; 
	tarSeq    = "."; 
	qrySeq    = "."; 
	strand    = '.'; 
	type      = "."; 
	score     = 0;   
	mismap    = 1.0;
	isClear   = false;

	return;
}

VarUnit::VarUnit (const VarUnit& V) {

	target = V.target; query   = V.query; tarSeq = V.tarSeq; 
	qrySeq = V.qrySeq; strand  = V.strand;
	type   = V.type  ; isClear = V.isClear; 

	score  = V.score ;  mismap = V.mismap;
	exp_target = V.exp_target;
	exp_tarSeq = V.exp_tarSeq;

	return;
}

void VarUnit::ConvQryCoordinate ( unsigned int qrySeqLen ) {

	// This funtion just conversion the coordinate of Axt/MAF format creat 
	// by 'lastz'/'last ', which mapped to the '-' strand
	if ( strand != '-' ) return;
	unsigned int itemp = query.start;
	query.start = qrySeqLen - query.end + 1;
	query.end   = qrySeqLen - itemp + 1;

	return;
}

void VarUnit::Swap(){ // Swap the target and query region. Ignore 'exp_target' 

	Region tmp = target; target = query;  query  = tmp;
	string str = tarSeq; tarSeq = qrySeq; qrySeq = str;

	return;
}

VarUnit VarUnit::ReAlign(Fa &targetSeq, Fa &querySeq, AgeOption opt){
// Return new VarUnit after AGE Realignment
// This is a design strategy, I'm not going to simply update the raw VarUnit!

	targetSeq.CheckFaId(target.id); // make sure target fa is right
	querySeq.CheckFaId(query.id);   // make sure query  fa is right

	AgeAlignment alignment(*(this), opt);
	if (alignment.Align(targetSeq.fa[target.id], querySeq.fa[query.id])){
	// Successful align!
		return alignment.vu();
	} else {
		return *(this);
	}
}

void VarUnit::OutStd(unsigned int tarSeqLen, unsigned int qrySeqLen, ofstream &O){ 
// Output the alignment to STDERR

	if (tarSeq.empty() || qrySeq.empty()){ 
		std::cerr << "tarSeq.empty() || qrySeq.empty()" << endl; exit(1); 
	}

	unsigned int qnl = NLength ( qrySeq );
	unsigned int tnl = NLength ( tarSeq );
	//O << target.id << "\t" << target.start << "\t" << target.end << "\t" 
	cout << target.id << "\t" << target.start << "\t" << target.end << "\t" 
	  << target.end - target.start + 1     << "\t" << double(tnl)/tarSeq.length()
	  << "\t" << tarSeqLen                 << "\t"
	  << query.id  << "\t" << query.start  << "\t" << query.end << "\t" 
	  << query.end  - query.start  + 1     << "\t" << double(qnl)/qrySeq.length()
	  << "\t" << qrySeqLen <<"\t"<< strand << "\t" << score << "\t" << mismap    
	  << "\t" << type      << endl;
	return;
}

void VarUnit::OutStd(unsigned int tarSeqLen, unsigned int exp_tarSeqLen, 
					 unsigned int qrySeqLen, ofstream &O){

	if (exp_target.isEmpty()) cerr << "[ERROR]exp_target is empty!\n";
	OutStd(tarSeqLen, qrySeqLen, O);

	if (exp_tarSeq.empty()){ cerr << "exp_tarSeq.empty() \n"; exit(1); }

	unsigned int qnl = NLength ( qrySeq     );
	unsigned int tnl = NLength ( exp_tarSeq );
	//O << exp_target.id << "\t" << exp_target.start     << "\t" << exp_target.end 
	cout << exp_target.id << "\t" << exp_target.start     << "\t" << exp_target.end 
	  << "\t" << exp_target.end - exp_target.start + 1 << "\t" 
	  << double(tnl)/exp_tarSeq.length()      << "\t" << exp_tarSeqLen << "\t" 
	  << query.id  << "\t"   << query.start   << "\t" << query.end     << "\t" 
	  << query.end  - query.start  + 1        << "\t" << double(qnl)/qrySeq.length() 
	  << "\t" << qrySeqLen  << "\t" << strand << "\t" << score << "\t" << mismap 
	  << "\t" << type + "-E"<< endl;
	return;
}

/**********************************
 * Author : Shujia Huang
 * Date   : 2014-08-05 22:49:56
 *
 * Class  AgeAlignment
 **********************************/

void AgeAlignment::Init(VarUnit &v, AgeOption opt){

	vu_     = v;
	para_   = opt;
	isInit_ = true;
}

void AgeAlignment::ExtendVariant(unsigned long int tarFaSize, 
								 unsigned long int qryFaSize,
								 int extandFlankSzie){
	
	if (!isInit_) { 
		cerr<<"[ERROR]You should init AgeAlignment before calling ExtendVariant()\n";
		exit(1);
	}
	
	vu_.target.start -= extandFlankSzie;
	vu_.target.end   += extandFlankSzie;
	vu_.query.start  -= extandFlankSzie;
	vu_.query.end    += extandFlankSzie;

	if (vu_.target.start < 1) vu_.target.start = 1;
    if (vu_.target.end   > tarFaSize) vu_.target.end = tarFaSize;	
    if (vu_.query.start  < 1) vu_.query.start = 1; 
    if (vu_.query.end    > qryFaSize) vu_.query.end  = qryFaSize;

	return;
}

bool AgeAlignment::IsHugeMemory(unsigned long int n, unsigned long int m){ 

	return (5 * n * m / 1000000000 > 10); // 10G
}

bool AgeAlignment::Align(string &tarFa, string &qryFa){

	if (!isInit_) { 
		cerr << "[ERROR] You should init AgeAlignment before calling Align()\n";
		exit(1);
	}

	int flag = 0;
    if (para_.indel) flag |= AGEaligner::INDEL_FLAG;
    if (para_.inv  ) flag |= AGEaligner::INVERSION_FLAG;
    if (para_.invl ) flag |= AGEaligner::INVL_FLAG;
    if (para_.invr ) flag |= AGEaligner::INVR_FLAG;
    if (para_.tdup ) flag |= AGEaligner::TDUPLICATION_FLAG;
    if (flag == 0) flag  = AGEaligner::INDEL_FLAG; // Default

    int n_modes = 0;
    if (flag & AGEaligner::INDEL_FLAG)        n_modes++;
    if (flag & AGEaligner::TDUPLICATION_FLAG) n_modes++;
    if (flag & AGEaligner::INVR_FLAG)         n_modes++;
    if (flag & AGEaligner::INVL_FLAG)         n_modes++;
    if (flag & AGEaligner::INVERSION_FLAG)    n_modes++;

    if (n_modes != 1) {
        cerr << "# [Error] In mode specification. ";
        if ( n_modes == 0 ) cerr << "No mode is specified.\n";
        if ( n_modes  > 1 ) cerr << "More than one mode is specified.\n";
        exit(1);
    }
	//////////////////////////////////////////////////////////////////////
	
	ExtendVariant(tarFa.length(), qryFa.length(), para_.extendVarFlankSzie);
    // Do not AGE if the memory cost is bigger than 10G.
    unsigned long int varTarSize = vu_.target.end - vu_.target.start;
	unsigned long int varQrySize = vu_.query.end  - vu_.query.start;
    if (IsHugeMemory(varTarSize, varQrySize)) return false;

#ifdef AGE_TIME
timeval ali_s,ali_e;
gettimeofday(&ali_s, NULL);
#endif

	Sequence *tarSeq = new Sequence(tarFa, vu_.target.id, 1, false);
	Sequence *qrySeq = new Sequence(qryFa, vu_.query.id , 1, false);
	Sequence *tar = tarSeq->substr(vu_.target.start, vu_.target.end);
	Sequence *qry = qrySeq->substr(vu_.query.start , vu_.query.end );
	Scorer scr(para_.match, para_.mismatch, para_.gapOpen, para_.gapExtend);

	bool isalign(true);
	if (para_.both) {
            
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
			//aligner1.printAlignment();
aligner1.Test(); 
        } else {
			//aligner2.printAlignment();
aligner2.Test();
        } 
        delete qryClone;

    } else {

        AGEaligner aligner(*tar, *qry);
        if (aligner.align(scr, flag)){
            //aligner.printAlignment();
aligner.Test();
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


