#include "VarUnit.h"

VarUnit::VarUnit() {

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

VarUnit::VarUnit(const VarUnit& V) {

	target = V.target; query   = V.query; tarSeq = V.tarSeq; 
	qrySeq = V.qrySeq; strand  = V.strand;
	type   = V.type  ; isClear = V.isClear; 

	score  = V.score ;  mismap = V.mismap;
	exp_target = V.exp_target;
	exp_tarSeq = V.exp_tarSeq;

	return;
}

void VarUnit::ConvQryCoordinate(unsigned int qrySeqLen) {

	// This funtion just conversion the coordinate of Axt/MAF format creat 
	// by 'lastz'/'last ', which mapped to the '-' strand
	if ( strand != '-' ) return;
	unsigned int itemp = query.start;
	query.start = qrySeqLen - query.end + 1;
	query.end   = qrySeqLen - itemp + 1;

	return;
}

void VarUnit::Swap() { // Swap the target and query region. Ignore 'exp_target' 

	Region tmp = target; target = query;  query  = tmp;
	string str = tarSeq; tarSeq = qrySeq; qrySeq = str;

	return;
}

vector<VarUnit> VarUnit::ReAlignAndReCallVar(Fa &targetSeq, Fa &querySeq, AgeOption opt) {
// Return new VarUnit after AGE Realignment
// This is a design strategy, I'm not going to simply update the raw VarUnit!

	targetSeq.CheckFaId(target.id); // make sure target fa is right
	querySeq.CheckFaId(query.id);   // make sure query  fa is right

	vector<VarUnit> vus;
	AgeAlignment alignment(*(this), opt);
	if (alignment.Align(targetSeq.fa[target.id], querySeq.fa[query.id])){
	// Successful align!
		vus = alignment.VarReCall();
	} 
	return vus;
}

void VarUnit::OutErr() {
// Output the alignment to STDERR

    if (tarSeq.empty() || qrySeq.empty()){
        std::cerr << "tarSeq.empty() || qrySeq.empty()" << endl; exit(1);
    }

    unsigned int qnl = NLength ( qrySeq );
    unsigned int tnl = NLength ( tarSeq );
    cerr << target.id << "\t" << target.start << "\t" << target.end << "\t"
      << target.end - target.start + 1     << "\t" << double(tnl)/tarSeq.length()
      << "\t" << query.id << "\t" << query.start  << "\t" << query.end << "\t"
      << query.end  - query.start  + 1     << "\t" << double(qnl)/qrySeq.length()
      << "\t" << strand << "\t" << score << "\t" << mismap
      << "\t" << type      << endl;
    return;
}

void VarUnit::OutStd(unsigned int tarSeqLen, unsigned int qrySeqLen, ofstream &O) { 
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
					 unsigned int qrySeqLen, ofstream &O) {

	if (exp_target.isEmpty()) cerr << "[ERROR]exp_target is empty!\n";
	OutStd(tarSeqLen, qrySeqLen, O);

	if (exp_tarSeq.empty()){ cerr << "exp_tarSeq.empty() \n"; exit(1); }

	unsigned int qnl = NLength ( qrySeq     );
	unsigned int tnl = NLength ( exp_tarSeq );
	//O << exp_target.id << "\t" << exp_target.start     << "\t" << exp_target.end 
	cout << exp_target.id << "\t" << exp_target.start  << "\t" << exp_target.end 
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

void AgeAlignment::Init(VarUnit &v, AgeOption opt) {

	vu_     = v;
	para_   = opt;
	isInit_ = true;
}

bool AgeAlignment::Align(string &tarFa, string &qryFa) {

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
	
	ExtendVU(tarFa.length(), qryFa.length(), para_.extendVarFlankSzie);
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

	isalign_ = true;
	if (para_.both) {
            
        Sequence *qryClone = qry->clone(); 
        qryClone->revcom();
        AGEaligner aligner1(*tar, *qry);
        AGEaligner aligner2(*tar, *qryClone);                                         
        bool res1 = aligner1.align(scr, flag);
		bool res2 = aligner2.align(scr, flag);

		if (!res1 && !res2) {
			cerr<<"No alignment made.\n";
			isalign_ = false;
		} else if (aligner1.score() >= aligner2.score()) {
			aligner1.SetAlignResult(); 
			alignResult_ = aligner1.align_result();
        } else {
			aligner2.SetAlignResult();
			alignResult_ = aligner2.align_result();
        } 
        delete qryClone;

    } else {

        AGEaligner aligner(*tar, *qry);
        if (aligner.align(scr, flag)){
			aligner.SetAlignResult();
			alignResult_ = aligner.align_result();
        } else {
            cerr << "No alignment made.\n";
			isalign_ = false; 
        }
    }

/*
if (isalign_) {
for (size_t i(0); i < alignResult_._map.size(); ++i) {
	cout << alignResult_._map[i].first._id << " " << alignResult_._map[i].first._start << "\t" << alignResult_._map[i].first._end << "\t" << alignResult_._map[i].first._sequence << "\n";
	cout << "Mapinfo:\t" << alignResult_._map_info[i] << "\n";
	cout << alignResult_._map[i].second._id << " " << alignResult_._map[i].second._start << "\t" << alignResult_._map[i].second._end << "\t" << alignResult_._map[i].second._sequence << "\n";
}
}
*/

#ifdef AGE_TIME
gettimeofday(&ali_e, NULL);
cout << "\nAlignment time is " << ali_e.tv_sec - ali_s.tv_sec
    + (ali_e.tv_usec - ali_s.tv_usec)/1e+6 <<" s\n\n";
#endif

	Sequence::deleteSequences(tarSeq);
	Sequence::deleteSequences(qrySeq);

	return isalign_;
}

vector<VarUnit> AgeAlignment::VarReCall() {
// debug here carefully

	vector<VarUnit> vus;
	if (isalign_) {

		isgoodAlign_ = false;
		if (alignResult_._identity.size() >= 3 &&
			alignResult_._identity[0].second > 98 && // Average identity
			alignResult_._identity[1].first  > 30 && // left side size
			alignResult_._identity[1].second > 95 && // left identity
			alignResult_._identity[2].first  > 30 && // right side size
			alignResult_._identity[2].second > 95)   // right identity 
			isgoodAlign_ = true;

		if (alignResult_._identity.size() > 3) {
			cerr << "[WARNING] alignResult_._identity.size() > 3!! It may be";
			cerr << "a bug, please contact the author to confirm!\n";
		}
			
		if (alignResult_._map.size() > 1) {
		// Call the variant in the excise region
			pair<MapData, MapData> pre_map = alignResult_._map[0];
			for (size_t i(1); i < alignResult_._map.size(); ++i) {
			// Variant in excise region	
				VarUnit var = CallVarInExcise(pre_map, alignResult_._map[i], 
											  alignResult_._strand);
				pre_map = alignResult_._map[i];
				vus.push_back(var);
			}
		}
		// Call the variant in the flank sequence of variant
		for (size_t i(0); i < alignResult_._map.size(); ++i) {
			vector<VarUnit> var = CallVarInFlank(alignResult_._map[i], 
												 alignResult_._map_info[i],
												 alignResult_._strand);
			for (size_t i(0); i < var.size(); ++i) vus.push_back(var[i]); 
		}
	} else {
		vus.push_back(vu_); // No this is not good for the "No alignment made" situation!!!
	}	

	return vus;	
}

VarUnit AgeAlignment::CallVarInExcise(pair<MapData, MapData> &lf, // Left side
									  pair<MapData, MapData> &rt, // Right side
									  char strand) {

	if (lf.first._id != rt.first._id || lf.second._id != rt.second._id) {
		cerr << "[BUG] The id is not match!!\n" << lf.first._id << ", "
			 << rt.first._id  << "; " << lf.second._id          << ", " 
			 << rt.second._id << "\n";
		exit(1);
	}

	int qlen = abs(rt.second._start - lf.second._end - 1); // Query
	int tlen = rt.first._start  - lf.first._end  - 1;      // Reference
	if (tlen < 0) cerr << "[WARNING] It'll cause bug below!\n";

	VarUnit vu;
	vu.strand    = strand;
	vu.score     = vu_.score; // It's the LAST aligne score
	vu.mismap    = vu_.mismap;// It's the LAST mismap probability
	vu.target.id = lf.first._id;
	vu.query.id  = lf.second._id;
	if (tlen > 0 && qlen == 0) {
	// Pure-Deletion
		vu.type = "DEL";
		// I should treat the type to be 'SDEL' 
		// if they're agree with the original type.
		if (toupper(vu_.type) == "DEL") vu.type = "SDEL";

		// Reference position
		vu.target.start= lf.first._end;
		vu.target.end  = vu.target.start + tlen;

		// Query position 
		// CAUTION : Just use the reference base at ALT field in VCF file
		vu.query.start = (strand == '+') ? lf.second._end : rt.second._start;
		vu.query.end   = vu.query.start;

	} else if (qlen > 0 && tlen == 0) {
	// Pure-Insertion
		vu.type = "INS";
		// I should treat the type to be 'SINS' 
        // if they're agree with the original type.
		if (toupper(vu_.type) == "INS") vu.type = "SINS";

		// Reference position
		vu.target.start = lf.first._end;
		vu.target.end   = vu.target.start;
        // Query position
		// CAUTION: Here is just the variant region, not include the position
		// which at the boundary of variant. So that I'll add one base of ref-
		// erence at ALT field in VCF file.
		vu.query.start = (strand == '+') ? lf.second._end + 1 : rt.second._start + 1;
		vu.query.end   = vu.query.start + qlen - 1;
	} else if (tlen == qlen) {
		vu.type = (tlen == 1) ? "SNP" : "BSubstitution";
		vu.target.start = lf.first._end + 1;
		vu.target.end   = vu.target.start + tlen -1;
		vu.query.start  = (strand == '+') ? lf.second._end + 1 : rt.second._start + 1;
		vu.query.end    = vu.query.start + qlen - 1;
	} else {
		// Simultaneous gap or Unknown Type
		// Actrually, "Unknown" should be imposible!!
		vu.type = (tlen != qlen && tlen > 0) ? "Sgap" : "Unknown" ;
		vu.target.start = lf.first._end;
		vu.target.end   = vu.target.start + tlen;
		vu.query.start  = (strand == '+') ? lf.second._end : rt.second._start;
		vu.query.end    = vu.query.start + qlen;
	}

vu.OutErr(); // Debug

	return vu;
}

vector<VarUnit> AgeAlignment::CallVarInFlank(pair<MapData, MapData> &m, 
                                             string &mapInfo, char strand) {


/**
 * GNAGGAGGTAGGCAGATCC-TGGGGCCAGTGGCATATGGGGCCTGGACACAGGGCGGCCT first
 * |.||||||||||||||||| |||||||||||||| |||||||||||||.||||||||||| Map Info
 * GGAGGAGGTAGGCAGATCCCTGGGGCCAGTGGCA-ATGGGGCCTGGACTCAGGGCGGCCT second
 **/

// Need Debug here!!!

	int inc1 = 1;
    int inc2 = 1; if (strand == '-') inc2 = -1;
	int pos1start(m.first._start), pos2start(m.second._start);
	int pos1end(m.first._start),   pos2end(m.second._start);

	vector<VarUnit> vus;
	VarUnit vuTmp;
	vuTmp.strand    = strand;
    vuTmp.score     = vu_.score; // It's the LAST aligne score
    vuTmp.mismap    = vu_.mismap;// It's the LAST mismap probability
    vuTmp.target.id = m.first._id;
    vuTmp.query.id  = m.second._id;

	for (size_t i(0); i < mapInfo.size(); ++i) {

		if (mapInfo[i] == '|') {
		// Homo block

			vuTmp.type   = "Homo";
			vuTmp.tarSeq = ".";
			vuTmp.qrySeq = ".";
			while (mapInfo[i] == '|' && i < mapInfo.size()) {
				pos1end += inc1;
				pos2end += inc2;
				++i;
			}
			--i; // Rock back 1 position
		} else if (mapInfo[i] == ' ') {
		// New Indel!
			if (m.first._sequence[i]  == '-') pos1start -= 1;
			if (m.second._sequence[i] == '-') pos2start -= 1;

			while (mapInfo[i] == ' ' && i < mapInfo.size()) {
				if (m.first._sequence[i]  != '-') pos1end += inc1;
				if (m.second._sequence[i] != '-') pos2end += inc2;
				++i;
			}
			--i; // Rock back 1 position
		} else if (mapInfo[i] == '.') {
		// New SNP or Map 'N'!
			if (toupper(m.first._sequence[i] == 'N')) {
				vuTmp.type   = "N";
				vuTmp.tarSeq = "N";
            	vuTmp.qrySeq = ".";
			} else if (toupper(m.second._sequence[i]) == 'N') {
				vuTmp.type   = "N";
				vuTmp.tarSeq = ".";
            	vuTmp.qrySeq = "N";
			} else {
				vuTmp.type   = "SNP";
				vuTmp.tarSeq = m.first._sequence[i]; // char2str(m.first._sequence[i])
            	vuTmp.qrySeq = m.second._sequence[i];
			}
		} else {
		// Who knows...
			cerr << "[ERROR] What is it?!" << mapInfo[i] 
				 << " Must be some error in AGE aligne.\nDetail : \n";
			cerr << m.first._sequence << "\t" << m.first._id  << ":" 
				 << m.first._start    << "-"  << m.first._end << "\n"
				 << mapInfo << "\n"
				 << m.second._sequence << "\t" << m.second._id  << ":"
                 << m.second._start    << "-"  << m.second._end << "\n";
			exit(1);
		}

		vuTmp.target.start = pos1start;
        vuTmp.target.end   = pos1end;
        vuTmp.query.start  = pos2start;
        vuTmp.query.end    = pos2end;

		vus.push_back(vuTmp);

vuTmp.OutErr(); // Debug

        pos1start = pos1end + 1;
        pos2start = pos2end + 1;
	}

	return vus;
}

void AgeAlignment::ExtendVU(unsigned long int tarFaSize, 
							unsigned long int qryFaSize,
							int extandFlankSzie) {
	
	if (!isInit_) { 
		cerr<<"[ERROR]You should init AgeAlignment before calling ExtendVU()\n";
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

bool AgeAlignment::IsHugeMemory(unsigned long int n, unsigned long int m) { 

	return (5 * n * m / 1000000000 > 10); // 10G
}

