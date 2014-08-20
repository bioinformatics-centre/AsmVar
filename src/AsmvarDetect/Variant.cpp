/* 
 * Author : Shujia Huang
 * Date   : 2013-10-2
 *
 */ 
#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <utility>
#include "Variant.h"

using namespace std;

void Variant::CallnSeq() {

	assert(tarSeq.length() == qrySeq.length());
	VarUnit tmpnseq; 
    tmpnseq.target.id = target.id;
    tmpnseq.query.id  = query.id;
    tmpnseq.strand    = strand;
	tmpnseq.score     = score; 
	tmpnseq.mismap    = mismap;
    tmpnseq.type      = "N";
	
	int tarPos(target.start), qryPos(query.start);
	int tarI(0), qryI(0);
	vector< VarUnit > tmpVarVector;
	for ( int i(0); i < tarSeq.length(); ++i ) {
	// e.g. : tarSeq = "-ab-c-t--",
	//        qrySeq = "aa-ccdcdt",
		tarPos += tarI; qryPos += qryI;
		tarI = tarSeq[i] == '-' ? 0 : 1;
		qryI = qrySeq[i] == '-' ? 0 : 1;
		if (tarI == 0 || qryI  == 0) continue; // target or query is '-' 
		if (toupper(tarSeq[i]) != 'N' && toupper(qrySeq[i]) != 'N') continue;

		tmpnseq.target.start = tarPos; tmpnseq.target.end = tarPos;
		tmpnseq.query.start  = qryPos; tmpnseq.query.end  = qryPos;
		if (toupper(tarSeq[i]) == 'N' && toupper(qrySeq[i]) != 'N') {
			tmpnseq.tarSeq = "N"; 
			tmpnseq.qrySeq = ".";
		} else if (toupper(tarSeq[i]) != 'N' && toupper(qrySeq[i]) == 'N') {
			tmpnseq.tarSeq = "."; 
			tmpnseq.qrySeq = "N";
		} else { // tarSeq[i] and qrySeq[i] are 'N'
			tmpnseq.tarSeq = "N";
			tmpnseq.qrySeq = "N";
		}
		tmpVarVector.push_back(tmpnseq);
	}
	tmpVarVector = MergeVarUnit( tmpVarVector );
	for (size_t i(0); i < tmpVarVector.size(); ++i) {
		// coordinates uniform to the positive strand
		tmpVarVector[i].ConvQryCoordinate(qryfa.fa[query.id].length());
		nSeq.push_back(tmpVarVector[i]);
	}

	return;
}

void Variant::CallHomoRef () {
    
    assert(tarSeq.length() == qrySeq.length());                                       
    VarUnit tmphseq; 
    tmphseq.target.id = target.id; 
    tmphseq.query.id  = query.id;
    tmphseq.strand    = strand;
    tmphseq.score     = score;
	tmphseq.mismap    = mismap; 
    tmphseq.type      = "Homo";

    int tarPos(target.start), qryPos(query.start); 
    int tarI(0), qryI(0);
    vector< VarUnit > tmpVarVector;
    for ( int i(0); i < tarSeq.length(); ++i ) {
    // e.g. : tarSeq = "-ab-c-t--", 
    //        qrySeq = "aa-ccdcdt",
        tarPos += tarI; qryPos += qryI;
        tarI = tarSeq[i] == '-' ? 0 : 1;
        qryI = qrySeq[i] == '-' ? 0 : 1;
        if (tarI == 0 || qryI  == 0) continue; // target or query is '-' 
		if (toupper(tarSeq[i]) == 'N' || toupper(qrySeq[i]) == 'N') continue;

		if ( toupper(tarSeq[i]) == toupper(qrySeq[i]) ) {
			tmphseq.target.start = tarPos; tmphseq.target.end = tarPos;
			tmphseq.query.start  = qryPos; tmphseq.query.end  = qryPos;
			tmphseq.tarSeq = ".";
			tmphseq.tarSeq = ".";
        	tmpVarVector.push_back(tmphseq);
		}
    }
    tmpVarVector = MergeVarUnit( tmpVarVector );
    for (size_t i(0); i < tmpVarVector.size(); ++i){
		// coordinates uniform to the positive strand!
		tmpVarVector[i].ConvQryCoordinate(qryfa.fa[query.id].length());
        homoRef.push_back(tmpVarVector[i]);
	}

    return;
}

void Variant::CallSNP () {

	assert ( tarSeq.length() == qrySeq.length() );

	VarUnit tmpsnp;
    tmpsnp.target.id = target.id;
    tmpsnp.query.id  = query.id;
    tmpsnp.strand    = strand;
	tmpsnp.score     = score;   // 2014-06-20 09:58:32
	tmpsnp.mismap    = mismap;  // 2014-06-20 09:58:39
    tmpsnp.type      = "SNP";

	int tarPos(target.start), qryPos(query.start);
	int tarI(0), qryI(0);
	for ( int i(0); i < tarSeq.length(); ++i ) {
	// e.g. : tarSeq = "-ab-c-t--",
	//        qrySeq = "aa-ccdcdt",
		tarPos += tarI;
		qryPos += qryI;
		tarI = tarSeq[i] == '-' ? 0 : 1;
		qryI = qrySeq[i] == '-' ? 0 : 1;
		if (tarI == 0 || qryI  == 0) continue; // target or query is '-'
		if (toupper(tarSeq[i]) == 'N' || toupper(qrySeq[i]) == 'N') continue;

		if ( toupper(tarSeq[i]) != toupper(qrySeq[i]) ) {
			tmpsnp.target.start = tarPos; tmpsnp.target.end = tarPos;
			tmpsnp.query.start  = qryPos; tmpsnp.query.end  = qryPos;
			tmpsnp.ConvQryCoordinate( qryfa.fa[query.id].length() ); // coordinates uniform to the positive strand!
			snp.push_back( tmpsnp );
		}
	}
}

void Variant::CallInsertion () { 
// Actually, we just have to call  the target gap regions.	
// All the coordinate of query should be uniform to the positive strand!
	vector< VarUnit > gap = CallGap ( target, tarSeq, query, qrySeq, strand, score, mismap, "INS" );
	for ( size_t i(0); i < gap.size(); ++i ) { 

		if ( !qryfa.fa.count(query.id) ) err("Missing some query id or query id can't match!!!\nThe unmatch query(main): " + query.id);

		gap[i].ConvQryCoordinate( qryfa.fa[query.id].length() ); // coordinates uniform to the positive strand!
		insertion.push_back( gap[i] );
		summary["query-Intra-gap('N')"] += qryfa.Nlength( query.id, gap[i].query.start, gap[i].query.end );
	}
}

void Variant::CallDeletion () { 
// Actually, we just have to call  the query gap regions.
// All the coordinate of query should be uniform to the positive strand!
	vector< VarUnit > gap = CallGap ( query, qrySeq, target, tarSeq, strand, score, mismap, "DEL" );
	for ( size_t i(0); i < gap.size(); ++i ) { 
		gap[i].Swap(); // Swap query and target region!
		if ( !qryfa.fa.count( query.id ) ) { err ("Missing some query id or query id can't match!!!\nThe unmatch query(main): "+query.id); } 
		gap[i].ConvQryCoordinate( qryfa.fa[query.id].length() ); // coordinates uniform to the positive strand!

		deletion.push_back( gap[i] );
		summary["query-Intra-gap('N')"] += qryfa.Nlength( query.id, gap[i].query.start, gap[i].query.end );
	}
}

bool Variant::CallIversion( MapReg left, MapReg middle, MapReg right ) {

	assert ( left.query.id == right.query.id && left.query.id == middle.query.id );

	bool flag(false);
	if ( left.target.id != right.target.id || left.target.id != middle.target.id ) return flag;

	unsigned int rStart = min ( left.target.start, right.target.start );
	unsigned int rEnd   = max ( left.target.end  , left.target.end    ); 
	if ( left.strand   != right.strand ||
		 middle.strand == left.strand  ||
		 middle.target.start < rStart  ||
		 middle.target.start > rEnd  ) return flag; // Just make sure the middle.target.start between [rStart,rEnd]

	VarUnit reg;
	flag       = true;
    reg.type   = "Inv";
    reg.target = middle.target;
	reg.query  = middle.query ;
    reg.strand = middle.strand;
	reg.score  = middle.score;  // 2014-06-20 10:35:21
    reg.mismap = middle.mismap; // 2014-06-20 10:35:15
	inversion.push_back( reg );

	return flag;
}

bool Variant::CallTranslocat ( MapReg left, MapReg middle, MapReg right ) {

	assert ( left.query.id == right.query.id && left.query.id  == middle.query.id );
	bool flag(false);
	// The target id should be the same of 'left' and 'right'
	if ( left.target.id != right.target.id ) return flag;

	unsigned int rStart = min ( left.target.start, right.target.start );
	unsigned int rEnd   = max ( left.target.end  , left.target.end    ); 

	// Do not overlap with 'left' or 'right'. And be the same strand on either side of the 'middle' region
	if ( left.strand != right.strand || (middle.target.start <= rEnd && middle.target.end >= rStart) ) return flag;

	flag = true;
	VarUnit reg, gap;

	gap = CallGap ( left, right );
	reg.exp_target = gap.target;

	reg.target = middle.target;
	reg.query  = middle.query ;
	reg.strand = middle.strand;
	reg.score  = middle.score;  // 2014-06-20 10:35:21
	reg.mismap = middle.mismap; // 2014-06-20 10:35:15
	reg.type   = ( middle.target.id == left.target.id ) ? "Trans-Intra" : "Trans-Inter";
	if ( middle.strand == left.strand ) {
		reg.type += "1";
	} else {
		reg.type += "2";
	}
	translocation.push_back( reg );

	return flag;
}

void Variant::GetMapReg () {

	// Call the query coverting function here to make the '-' strand coordinates 
	// of query be the same as '+' strand.
	// ConvQryCoordinate(...) is a memerber function of class 'MAF'  
	qryfa.CheckFaId(query.id);
	ConvQryCoordinate(qryfa.fa[query.id].length());

	MapReg mpI;
	mpI.target = target; mpI.query = query; mpI.strand = strand;
	mpI.score  = score ; mpI.mismap= mismap; // 2014-06-20 10:39:54
	mapreg[query.id].push_back( mpI );       // The 'key' is query.id

	target.info = query.id;
	query.info  = target.id;
	maptar[target.id].push_back( target ); // stored the mapped target regions here
	mapqry[query.id].push_back ( query  ); // stored the mapped query  regions here
}

void Variant::Assign2allvariant(vector<VarUnit> &v) {

	for (size_t i(0); i < v.size(); ++i) {
		allvariant[v[i].target.id].push_back(v[i]);
	}
	return;
}

void Variant::AGE_Realign() {

	allvariant.clear(); // make sure the value is empty

	//AGE_Realign(snp); // No need!!
	//AGE_Realign(intragap); // No need!
	AGE_Realign(insertion);  // New variant will store in 'allvariant' 
	AGE_Realign(deletion);   // New variant will store in 'allvariant' 
	AGE_Realign(inversion);  // New variant will store in 'allvariant' 
	//AGE_Realign(translocation);
	AGE_Realign(simulreg);   // New variant will store in 'allvariant' 
	AGE_Realign(nosolution); // New variant will store in 'allvariant' 
	//AGE_Realign(clipreg); No need!
	//AGE_Realign(nomadic); No need!

	Assign2allvariant(snp);     // store in 'allvariant'
	Assign2allvariant(homoRef); // store in 'allvariant'
	Assign2allvariant(nSeq);    // store in 'allvariant'
	
	//Assign2allvariant(translocation);

	map<string, size_t> index;
	for (map<string, vector<VarUnit> >::iterator it(allvariant.begin()); 
		it != allvariant.end(); ++it) { 

		Unique(it->second); 
		sort(it->second.begin(), it->second.end(), MySortByTarV);
		index[it->first] = 0;
	}

	return;
}

void Variant::AGE_Realign(vector<VarUnit> &R) {

	// re-aligne :
	AgeOption opt;
	VarUnit vu;
	vector<VarUnit> vus;
	for (size_t i(0); i < R.size(); ++i) {

		// R[i] should be replace by 'v' after ReAlign!
		vector<VarUnit> v = R[i].ReAlignAndReCallVar(tarfa, qryfa, opt);
		// Not going to deal with the flankin region
		if (v.empty()) continue;
		if (v[0].type.find("-AGE") == string::npos) { // has variant in exci-reg
			if (R[i].Empty()) v[0].Clear(); // Can just happen after call Filter()
			allvariant[v[0].target.id].push_back(v[0]);  
		}
		//for (size_t j(0); j < v.size(); ++j) vus.push_back(v[j]);
		
cerr << "\n***********************************\n";
R[i].OutErr();
cerr << "\n********** AGE Process ************\n";
//for (size_t j(0); j < v.size(); ++j) v[j].OutErr();
v[0].OutErr();
cerr << "## allvariant ## " << itoa(allvariant[v[0].target.id].back().cipos.first) << "\n";
allvariant[v[0].target.id].back().OutErr();
cerr << "******* GOOD ******************\n";
	}

	return;
}

void Variant::Unique(vector<VarUnit> &v) {

	cerr << "#[INFO] Masking all the repeat appear varaints.\n";
	set<string> hasAppear;
	for (size_t i(0); i < v.size(); ++i) {
		string key = v[i].target.id + ":" + itoa(v[i].target.start) + ":"
					+ itoa(v[i].target.end) + "-" + v[i].query.id   + ":"
					+ itoa(v[i].query.start)+ ":" + itoa(v[i].query.end);
		if (hasAppear.count(key)) v[i].Clear(); // Mask the repeat appear var!
		hasAppear.insert(key);
	}
}

map< string, vector<Region> > Variant::VarTarRegs() {

	map< string, vector<Region> > varTarRegs;

	for (size_t i(0); i < insertion.size(); ++i ) {
		if ( insertion[i].Empty() ) continue;
		insertion[i].target.info = insertion[i].query.id;
		varTarRegs[insertion[i].target.id].push_back(insertion[i].target);
	}
	for (size_t i(0); i < deletion.size(); ++i ) {
        if ( deletion[i].Empty() ) continue;
		deletion[i].target.info = deletion[i].query.id;
        varTarRegs[deletion[i].target.id].push_back(deletion[i].target);
    }
	for (size_t i(0); i < simulreg.size(); ++i ) {
		simulreg[i].target.info = simulreg[i].query.id;
		varTarRegs[simulreg[i].target.id].push_back(simulreg[i].target);
	}

	for ( map< string,vector<Region> >::iterator it( varTarRegs.begin() ); it != varTarRegs.end(); ++it ) {
        sort ( varTarRegs[it->first].begin(), varTarRegs[it->first].end(), SortRegion );
    }

	return varTarRegs;	
}

void Variant::CallSV () { 
// Call Stuctural variants, not indels!!! 
// Just use the memerber value 'mapreg' in this memerber function.
// Should be debug carfully here!
// All the coordinate of query should be uniform to the positive strand, then I can sort them!

	map< size_t, vector<size_t> > mark;
	VarUnit simulgap;
	for ( map<string, vector<MapReg> >::iterator it( mapreg.begin() ); it != mapreg.end(); ++it ) {
		if ( it->second.size() < 2 ) continue; // More than 2

		map<unsigned int, size_t> qry2tarOrder;
		sort( it->second.begin(), it->second.end(), MySortByTarM ); // Sort by the coordinate of target mapping positions
		for ( size_t i(0); i < it->second.size(); ++i ) { qry2tarOrder[it->second[i].query.start] = i; }
		sort( it->second.begin(), it->second.end(), MySortByQryM ); // Sort by the coordinate of query  mapping positions

		vector<MapReg> tmpreg;
		unsigned int pos1, pos2, pos;
		for ( size_t i(0); i < it->second.size(); ++i ) {

			pos = it->second[i].query.start;   // right (side)
			if ( tmpreg.size() == 1 ) {
				pos1 = tmpreg[0].query.start;
				if ( labs(qry2tarOrder[pos] - qry2tarOrder[pos1]) == 1 ) {
					// Regular simultaneous gap regions
					if ( CallSimultan(tmpreg[0], it->second[i]) ) tmpreg.clear();
				}
			} else if ( tmpreg.size() == 2 ) { // Candidate translocation or Candidate Inversion. It's not success when calling simultaneous gap
				pos1 = tmpreg[0].query.start;  // left 
				pos2 = tmpreg[1].query.start;  // middle
				if ( labs(qry2tarOrder[pos] - qry2tarOrder[pos1]) == 1 ) { // Means tmpreg[0] and tmpreg[1] are not the neighbour
					
					if ( CallTranslocat(tmpreg[0], tmpreg[1], it->second[i]) ) { 
						tmpreg.clear(); // Here 'tmpreg' just contain the two head element, in fact, I'm just clear the first and second elements.
					} else {
						tmpreg.erase( tmpreg.begin() );
						if ( labs(qry2tarOrder[pos] - qry2tarOrder[pos2]) == 1 ) { 
							if ( CallSimultan(tmpreg[0], it->second[i]) ) tmpreg.clear(); // Call simultaneous gap
						}
					}
				} else if ( (labs(qry2tarOrder[pos2] - qry2tarOrder[pos1]) == 1) && (labs(qry2tarOrder[pos] - qry2tarOrder[pos2]) == 1) ){ 

					if ( CallIversion(tmpreg[0], tmpreg[1], it->second[i]) ) { 
						tmpreg.clear(); // Here 'tmpreg' just contain the two head element , in fact, I'm just clear the first and second elements.
					} else {
						tmpreg.erase ( tmpreg.begin() );
						if ( CallSimultan(tmpreg[0], it->second[i]) ) tmpreg.clear();
					}
				} else if (labs(qry2tarOrder[pos2] - qry2tarOrder[pos1]) == 1) { // No solution
					CallReg ( tmpreg[0], "Nos-JustNo", nosolution );
					tmpreg.erase ( tmpreg.begin() );
				} else {
					CallReg ( tmpreg[0], "Nos-Complex", nosolution );
					tmpreg.erase ( tmpreg.begin() );
				}
			} else if ( tmpreg.size() > 2 ) {
				std::cerr << "\n[ERROR] Program bugs! " << endl;
			}
			tmpreg.push_back( it->second[i] );
		}
		if ( tmpreg.size() == 2 ) {
			if ( CallSimultan ( tmpreg[0], tmpreg[1] ) ) {
				tmpreg.clear();
			} else if ( tmpreg[0].target.id == tmpreg[1].target.id ) { 
				CallReg ( tmpreg[0], "Nos-Strand", nosolution );
                CallReg ( tmpreg[1], "Nos-Strand", nosolution );
			} else { // Not in the same target
				CallReg ( tmpreg[0], "Nos-Complex", nosolution );
				CallReg ( tmpreg[1], "Nos-Complex", nosolution );
			}
		}
	}
}

bool Variant::CallSimultan ( MapReg left, MapReg right ) {

	bool success ( false );
	VarUnit simulgap;
	if ( left.target.id != right.target.id || left.strand != right.strand ) return success;	// (2013-10-20 19:17:40)

	simulgap = CallGap( left, right );
	simulreg.push_back( simulgap );
	success  = true;

	return success;
}

void Variant::CallReg ( MapReg mapreg, string type, vector< VarUnit > & varReg ) {

	VarUnit reg;
	reg.target = mapreg.target;
	reg.query  = mapreg.query;
	reg.strand = mapreg.strand;
	reg.score  = mapreg.score;  // 2014-06-20 10:37:45
	reg.mismap = mapreg.mismap; // 2014-06-20 10:37:55
	reg.type   = type;      // No solution
	varReg.push_back( reg );
}

void Variant::CallClipReg () {

	if ( qryfa.fa.empty() ) { std::cerr << "[WARNING] No Clip regions. Because the query fa is empty!" << endl; return; }
	VarUnit tmp;
	unsigned int rStart, rEnd;
	for ( map< string, vector<Region> >::iterator it( mapqry.begin() ); it != mapqry.end(); ++it ) {

		rStart = RegionMin( it->second );
		rEnd   = RegionMax( it->second );

		if ( rStart == 0 || rEnd == 0 ) { 

			std::cerr << "rStart == 0 || rEnd == 0" << "\nrStart: " << rStart << "\trEnd: " << rEnd << endl;
			for ( size_t i(0); i < it->second.size(); ++i )	
				std::cerr << "Size: " << it->second.size() << "\tit->second: " << it->second[i].id 
					 << "\t" << it->second[i].start   << "\t" << it->second[i].end << endl;
			exit(1);
		}

		if ( rStart == 1 && rEnd == qryfa.fa[it->first].length() ) ++summary["Full-Align"];

		tmp.query.id = it->first;
		tmp.strand   = '.';
		tmp.type     = "Clip";
		if ( rStart > 1 ) { tmp.query.start = 1; tmp.query.end = rStart - 1; clipreg.push_back( tmp ); }
		if ( rEnd   < qryfa.fa[it->first].length() ) {
			tmp.query.start = rEnd + 1;
			tmp.query.end   = qryfa.fa[it->first].length();
			clipreg.push_back( tmp );
		}
	}
	return;
}

void Variant::CallNomadic () {

	if (qryfa.fa.empty()) { std::cerr << "[WARNING] No Nomadic regions. Because the query fa is empty!" << endl; return; }
	VarUnit tmp;
	for ( map<string, string>::iterator it( qryfa.fa.begin() ); it != qryfa.fa.end(); ++it ) {

		if ( mapqry.count( it->first ) ) continue;
		tmp.query.id    = it->first;
		tmp.query.start = 1;
		tmp.query.end   = it->second.length();
		tmp.strand      = '.';
		tmp.type        = "Nomadic";
		nomadic.push_back( tmp );
	}
	return;
}

bool Variant::IsSameStrand ( vector<MapReg> & mapreg ) {

	bool same ( true );
	char strand = mapreg[0].strand;
	for ( size_t i(1); i < mapreg.size(); ++i ) {
		if ( strand != mapreg[i].strand ) { same = false; break; }
	}
	return same;
}

void Variant::Filter () {

	map< string,vector<Region> > filterReg;
	map<string, size_t> index;
	for ( size_t i(0); i < nosolution.size(); ++i ) 
		filterReg[nosolution[i].query.id].push_back( nosolution[i].query );

	for ( map< string,vector<Region> >::iterator it( filterReg.begin() ); 
		  it != filterReg.end(); ++it ){

		sort(filterReg[it->first].begin(), filterReg[it->first].end(), SortRegion);
		index[it->first] = 0;
	}

	// Filter insertion regions which in filterReg 
	FilterReg ( filterReg, index, insertion );
	// Filter deletion regions which in filterReg
	FilterReg ( filterReg, index, deletion  );
	// Filter SNP which in filterReg
	FilterReg ( filterReg, index, snp       );

	if ( snp.empty() ) return;
	map <string, vector<size_t> > tmp;
	for ( size_t i(0); i < snp.size(); ++i ) {
		string key = snp[i].target.id + ":" + snp[i].query.id + ":" + itoa( snp[i].target.start );
		tmp[key].push_back( i );
	}
	for ( map <string, vector<size_t> >::iterator it( tmp.begin() ); it != tmp.end(); ++it ) {
		if ( it->second.size() == 1 ) continue;
		for ( size_t i(0); i < it->second.size(); ++i ) snp[it->second[i]].Clear();
	}

	return;
}

void Variant::FilterReg(map< string,vector<Region> > tarregion, map<string, size_t> index, vector<VarUnit> &region) {
// Just used in Filter() function

	if ( region.empty() || tarregion.empty() ) return;
	sort ( region.begin(), region.end(), MySortByQryV );	

	string key;
	bool flag;
	unsigned int prepos = region[0].query.start; 
	for ( size_t i(0); i < region.size(); ++i ) {

		if ( !index.count(key) ) continue;

		assert( prepos <= region[i].query.start ); // Check sorted
		key = region[i].query.id;
		flag= true;

		for ( size_t j( index[key] ); j < tarregion[key].size(); ++j ) {
			if ( j > 0 ) assert( tarregion[key][j-1].start <= tarregion[key][j].start ); // Check sorted
			if ( region[i].query.start > tarregion[key][j].end   ) continue;
			if ( region[i].query.end   < tarregion[key][j].start ) break;

			if ( flag ) { flag = false; index[key] = j; }
			region[i].Clear();
		}
	}
	return;
}

void Variant::Summary( string file ) {

	summary["SNP"]     = snp.size();
	summary["Ins"]     = insertion.size();
	summary["Del"]     = deletion.size();
	summary["Inv"]     = inversion.size();
	summary["Clip"]    = clipreg.size();
	summary["Nomadic"] = nomadic.size();
	for ( size_t i(0); i < translocation.size(); ++i ) summary[translocation[i].type]++;
	for ( size_t i(0); i < nosolution.size();    ++i ) summary[nosolution[i].type]++;
	for ( size_t i(0); i < simulreg.size();      ++i ) summary[simulreg[i].type]++;

	map< string, unsigned int > tarCov, qryCov;
	map< string, vector<Region> >::iterator p( mapqry.begin() );
	for ( ; p != mapqry.end(); ++p ) { summary["qryCovlength"] += Covlength( p->second ); qryCov[p->first] += Covlength( p->second ); }
	p = maptar.begin();
	for ( ; p != maptar.end(); ++p ) { summary["tarCovlength"] += Covlength( p->second ); tarCov[p->first] += Covlength( p->second ); }

	ofstream O ( file.c_str() );
    if ( !O ) { std::cerr << "Cannot write to file : " << file << endl; exit(1); }
	O << "Summary Information for " << sample << "\n\n";
	for ( map<string, unsigned int>::iterator pt( summary.begin() ); pt!= summary.end(); ++pt ) 
		O << pt->first << "\t" << pt->second << "\n";
	
	O << "QryCovlength/querylength  " << double ( summary["qryCovlength"] ) / qryfa.length << "\n";
	O << "TarCovlength/targetlength " << double ( summary["tarCovlength"] ) / tarfa.length << "\n";
	O << "TarCovlength/targetlength(NO 'N') "<< double(summary["tarCovlength"])/(tarfa.length-tarfa.nsize) << "\n";
	O << "SNP/querylength           " << double(summary["SNP"]) / qryfa.length  << "\n";
	O << "SNP/targetlength          " << double(summary["SNP"]) / tarfa.length  << "\n";
	O << "\n";
	for ( map<string, unsigned int>::iterator p(tarCov.begin()); p != tarCov.end(); ++p ) {

		double ratio = double(p->second)/tarfa.fa[p->first].length();
		O << "# " << p->first << "\t" << p->second << "/" << tarfa.fa[p->first].length() << "\t" << ratio << "\n";
	}
	O << endl;
	O.close();
}

void Variant::Output ( string file ) {

	ofstream O ( file.c_str() );
    if ( !O ) { std::cerr << "Cannot write to file : " << file << endl; exit(1); }

	Output ( insertion, O ); //
	Output ( deletion,  O ); //
	Output ( inversion, O );
	Output ( simulreg,  O ); //
	Output ( clipreg,   O ); //
	Output ( nomadic,       O ); //
	Output ( nosolution,    O ); 
	Output ( translocation, O );

	Output ( snp,     O );
	Output ( homoRef, O );
	Output ( nSeq,    O );

	O.close();
	return;
}

void Variant::Output ( vector< VarUnit > & R, ofstream& O ) {

	sort (R.begin(), R.end(), MySortByTarV);
	for ( size_t i(0); i < R.size(); ++i ) {

		if ( R[i].Empty() ) continue;

		R[i].tarSeq = ( R[i].target.id == "-" ) ? "-" : tarfa.fa[R[i].target.id].substr(R[i].target.start - 1, R[i].target.end - R[i].target.start + 1);
		R[i].qrySeq = qryfa.fa[ R[i].query.id].substr( R[i].query.start - 1, R[i].query.end - R[i].query.start + 1 );
		if ( R[i].strand == '-' ) R[i].qrySeq = ReverseAndComplementary( R[i].qrySeq );

		if ( R[i].exp_target.isEmpty() ) {
			R[i].OutStd(  tarfa.fa[R[i].target.id].length(), qryfa.fa[ R[i].query.id].length(), O );
		} else {
			if ( R[i].exp_target.id.empty() || R[i].exp_target.id == "-" ) err ( "exp_target.id.empty() || R[i].exp_target.id == '-' " );
			R[i].exp_tarSeq = tarfa.fa[ R[i].exp_target.id ].substr( R[i].exp_target.start - 1, R[i].exp_target.end - R[i].exp_target.start + 1); 
			R[i].OutStd( tarfa.fa[R[i].target.id].length(), tarfa.fa[R[i].exp_target.id].length(), qryfa.fa[ R[i].query.id].length(), O);
		}
/*
// re-aligne :
AgeOption opt;
vector<VarUnit> vus = R[i].ReAlignAndReCallVar(tarfa, qryfa, opt);
cerr << "\n***********************************\n";
R[i].OutErr();
cerr << "\n********** AGE Process ************\n";
for (size_t i(0); i < vus.size(); ++i) vus[i].OutErr();
*/
	}
}

void Variant::OutputSNP ( string file ) {

	ofstream O ( file.c_str() );
	O <<  
"##fileformat=VCFv4.1                                                       \n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">              \n\
##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Mappind strand\">        \n\
##FORMAT=<ID=VT,Number=1,Type=String,Description=\"Variant type\">          \n\
##FORMAT=<ID=QR,Number=1,Type=String,Description=\"Query Info\">            \n\
##FORMAT=<ID=MS,Number=2,Type=Integer,Description=\"Mapping Score\">        \n\
##FORMAT=<ID=MIP,Number=1,Type=Float,Description=\"Mismapped probability\"> \n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n";
	sort (snp.begin(), snp.end(), MySortByTarV);
    for ( size_t i(0); i < snp.size(); ++i ) {
        if ( snp[i].Empty() ) continue;

        snp[i].tarSeq = ( snp[i].target.id == "-" ) ? "-" : tarfa.fa[snp[i].target.id].substr(snp[i].target.start - 1, snp[i].target.end - snp[i].target.start + 1);
// WSKMYRVDBHwskmyrvdbh
// ACTACAAATAactacaaata
		bool isNext (false);
		switch ( toupper(snp[i].tarSeq[0]) ) {
			case 'W' : ; case 'M':; case 'R':; case 'V':; case 'D':; 
			case 'H' : ; case 'S':; case 'Y':; case 'K':; 
			case 'B' : isNext = true; break;
		}

        snp[i].qrySeq = qryfa.fa[snp[i].query.id].substr( snp[i].query.start - 1, snp[i].query.end - snp[i].query.start + 1 );
		if ( snp[i].strand == '-' ) snp[i].qrySeq = ReverseAndComplementary( snp[i].qrySeq );
		if ( snp[i].tarSeq == snp[i].qrySeq ) continue;

		O << snp[i].target.id << "\t" << snp[i].target.start << "\t.\t" << snp[i].tarSeq << "\t" << snp[i].qrySeq << "\t255\tPASS\t.\tGT:ST:VT:QR:MS:MIP\t"
          << "./.:" << snp[i].strand  << ":" + snp[i].type + ":" + snp[i].query.id + "-" + itoa(snp[i].query.start) + "-" + itoa(snp[i].query.end) + ":"
          + itoa( snp[i].score ) + ":" << snp[i].mismap << "\n";
    }
	O.close();
}

void Variant::OutputGap( string file ) {
// The inter-gaps between different scaffolds in the same chromosome need to calculate.
// Abort : 2013-11-04 13:20:30 !!!

	ofstream O ( file.c_str() );
	if ( !O ) { std::cerr << "Cannot write to file : " << file << endl; exit(1); }

	map<string, vector<MapReg> > tmpmapreg;
	for ( map<string, vector<MapReg> >::iterator it( mapreg.begin() ); it != mapreg.end(); ++it ) {
		for ( size_t i(0); i < it->second.size(); ++i ) {
			tmpmapreg[it->second[i].target.id].push_back( it->second[i] );
		}
	}
	// Get inter scaffold gaps
	for ( map<string, vector<MapReg> >::iterator it( tmpmapreg.begin() ); it != tmpmapreg.end(); ++it ) {

		sort( it->second.begin(), it->second.end(), MySortByTarM ); // Sort by the coordinate of target mapping positions
		// Get Inter scaffold gaps' regions
		MapReg tmpMR = it->second[0];
		for ( size_t i(1); i < it->second.size(); ++i ) {
			//it->second[i].OutErrAlg();  // For Test
			if ( tmpMR.query.id == it->second[i].query.id ) {
				if (tmpMR.target.end < it->second[i].target.end) tmpMR = it->second[i];
				continue;
			}

			int length = (it->second[i].target.start - 1) - (tmpMR.target.end + 1) + 1;
			if (tmpMR.target.end >= it->second[i].target.start - 1) {
				O << tmpMR.target.id << "\t" << tmpMR.target.end  << "\t" << it->second[i].target.start << "\t" << length << "\tOverlap\t"
				  << tmpMR.target.id + ":"   << tmpMR.target.start<< "-"  << tmpMR.target.end << "|" 
				  << tmpMR.query.id  + ":"   << tmpMR.query.start << "-"  << tmpMR.query.end  << "#"
				  << it->second[i].target.id + ":" << it->second[i].target.start << "-" << it->second[i].target.end << "|"
				  << it->second[i].query.id  + ":" << it->second[i].query.start  << "-" << it->second[i].query.end  << "\n";
			} else {
				O << tmpMR.target.id << "\t" << tmpMR.target.end + 1 << "\t" << it->second[i].target.start - 1 << "\t" << length << "\tInterGap\t"
				  << tmpMR.target.id + ":"   << tmpMR.target.start   << "-"  << tmpMR.target.end << "|" 
                  << tmpMR.query.id  + ":"   << tmpMR.query.start    << "-"  << tmpMR.query.end  << "#"
				  << it->second[i].target.id + ":" << it->second[i].target.start << "-" << it->second[i].target.end << "|"
                  << it->second[i].query.id  + ":" << it->second[i].query.start  << "-" << it->second[i].query.end  << "\n";
				summary["InterScaffoldGap"] += length;
			}
			if (tmpMR.target.end < it->second[i].target.end) tmpMR = it->second[i];
		}
	}
	O.close();
}

// friendship function in class 'Variant'
unsigned int Variant::Covlength ( vector<Region> mapreg ) {

	if ( mapreg.empty() ) return 0;

	sort ( mapreg.begin(), mapreg.end(), SortRegion );
	unsigned int length(mapreg[0].end - mapreg[0].start + 1), prepos( mapreg[0].end );

	for ( size_t i(1); i < mapreg.size(); ++i ) {

		if ( mapreg[i].start > prepos ) {
			length += ( mapreg[i].end - mapreg[i].start + 1 );
			prepos = mapreg[i].end;
		} else if ( mapreg[i].end > prepos ) {
			length += ( mapreg[i].end - prepos );
			prepos = mapreg[i].end;
		} else {
			length += 0;
		}
	}
	
	return length;
}

vector< VarUnit > Variant::CallGap ( Region & tar,    // chrM 16308 16389  or Contig102837 1 81
                                     string & tarSeq, // CATAGTACATAAAGTCATTTACCGTACATAGCACATTACAG
                                     Region & qry,    // Contig102837 1 81 or chrM 16308 16389
                                     string & qrySeq, // CATAGTACATAAAGTCATTTACCGTACATAGCACATTACAG
                                     char   strand,   // + or -
                                     long   score,    // 221 
                                     double mismap,   // 1e-10
                                     string type ) {  // insertion or deletion

// This function is just used to call the gap regions of 'tar' (not for 'qry'!!). which will be indel actually ( indels are gaps ).
	assert ( tarSeq.length() == qrySeq.length() );

	VarUnit tmpgap; 
	tmpgap.target.id = tar.id;
	tmpgap.query.id  = qry.id;
	tmpgap.strand    = strand; 
	tmpgap.score     = score;
	tmpgap.mismap    = mismap;
	tmpgap.type      = type;

	vector< VarUnit > gap;
	size_t i(0), j(0), tgaplen(0), qgaplen(0);
	while ( ( i = tarSeq.find_first_of('-',j) ) != string::npos ) {
	// I do have a program to debug this part at : /ifs1/ST_EPI/USER/huangshujia/bin/cpp_bin/learn_cpp/test_string.cpp
	// e.g. : tarSeq = "-ab-c-t--", 
	//        qrySeq = "aa-ccdcdt",
		for ( size_t k(j); k < i; ++k ) { if( qrySeq[k] == '-' ) ++qgaplen; }
		tmpgap.query.start = qry.start + i - qgaplen; // Get the query start position which in the tar-gap region.

		if ( ( j = tarSeq.find_first_not_of('-',i) ) == string::npos ) j = tarSeq.length();
		for ( size_t k(i); k < j; ++k ) { if( qrySeq[k] == '-' ) ++qgaplen; }
		tmpgap.query.end  = qry.start + j - qgaplen - 1;      // Get the end  position which in the tar-gap region.

		tgaplen += ( j - i ); // Carefull, You must stop if catch 'j = tarSeq.length()' one time
		tmpgap.target.start = tar.start + j - tgaplen - 1;
		tmpgap.target.end   = tmpgap.target.start;
		gap.push_back( tmpgap );
	}

	return gap;
}

VarUnit Variant::CallGap ( MapReg left, MapReg right ) { 
// Call the simultaneous gap between 'left' mapped region and the 'right' one
// 'left' and 'right' should be the same target id, the same query id and the same strand!
	assert ( left.target.id == right.target.id && left.query.id == right.query.id && left.strand == right.strand );

	VarUnit gap;
	gap.target.id = left.target.id;
	gap.query.id  = left.query.id ;
	gap.strand    = left.strand;
	gap.score     = min(left.score , right.score ); // 2014-06-20 10:47:58
	gap.mismap    = min(left.mismap, right.mismap); // 2014-06-20 10:48:06

	gap.query.start = left.query.end;
	gap.query.end   = right.query.start;
	if ( left.target.start <= right.target.start ) { // Ascending

		// Should +1? what if the size larger than the size of query after +1. 
		// so I don't want to +1!
		gap.target.start = left.target.end;

		// I don't want to -1 ! If left.query overlap with right.query, we will
		// find gap.query.start >= gap.query.end, that's not right!
		gap.target.end   = right.target.start;
	} else { // Decending I consider the '-' strand here!!
		gap.target.start = right.target.end;
		gap.target.end   = left.target.start;
	}

	bool isTarOvlp( false ), isQryOvlp( false );
	if ( gap.target.end <= gap.target.start ) {
		gap.target.start = min ( left.target.start, right.target.start );
		gap.target.end   = max ( left.target.end  , right.target.end   );
		isTarOvlp        = true;
	}
	if ( gap.query.end  <= gap.query.start ) {
		gap.query.start  = min ( left.query.start, right.query.start );
		gap.query.end    = max ( left.query.end  , right.query.end   );
		isQryOvlp        = true;
	}

	if ( isTarOvlp && isQryOvlp ) { 
		gap.type = "gap-X1";
	} else if ( isQryOvlp ) {
		gap.type = "gap-X2";
	} else if ( isTarOvlp ) {
		gap.type = "gap-1";
	} else { // All not overlap
		gap.type = "gap-2"; // This gap means the normal simultaneous gaps
	}

	return gap;
}

void Variant::Output2VCF ( string file ) {

	VcfHeader header;
	string h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample;
	header.Add("#CHROM", h);

    ofstream O (file.c_str());
	map<string, string> hdata = header.data();
	for (map<string, string>::iterator it(hdata.begin()); 
		it != hdata.end(); ++it) {
		O << it->second << "\n";
	}
	
	for (map<int, string>::iterator it(tarfa.order2id.begin()); 
		 it != tarfa.order2id.end(); ++it) {

		if (!allvariant.count(it->second)) continue;// No such Fa id
		for (size_t i(0); i < allvariant[it->second].size(); ++i) {

			if (allvariant[it->second][i].type != "N" && 
				allvariant[it->second][i].type.find("Homo") == string::npos) {

				unsigned long int start = allvariant[it->second][i].target.start;
				unsigned long int end   = allvariant[it->second][i].target.end;
				allvariant[it->second][i].tarSeq = tarfa.fa[it->second].substr(
						start - 1, end - start + 1);

				start = allvariant[it->second][i].query.start;
				end   = allvariant[it->second][i].query.end;
				allvariant[it->second][i].qrySeq = 
					qryfa.fa[allvariant[it->second][i].query.id].substr(
						start - 1, end - start + 1);
			
				if (allvariant[it->second][i].strand == '-') 
					allvariant[it->second][i].qrySeq = ReverseAndComplementary(
											  allvariant[it->second][i].qrySeq);
			}

			VCF vcfline;
			vcfline.chrom_ = it->second;
			vcfline.pos_   = allvariant[it->second][i].target.start;
			vcfline.Id_    = ".";
			vcfline.ref_   = allvariant[it->second][i].tarSeq;
			vcfline.alt_   = allvariant[it->second][i].qrySeq;
			vcfline.qual_  = 255;
			VcfFormat format;
			if (allvariant[it->second][i].type == "N") {
				vcfline.filters_ = "NCALL";
				format.Set("GT", "./.");
			} else if (allvariant[it->second][i].type.find("Homo") != string::npos) {
				vcfline.filters_ = "REFCALL"; // Homo reference region
				format.Set("GT", "0/0");
			} else {
				vcfline.filters_ = (allvariant[it->second][i].isGoodReAlign) ? 
									"PASS-AGE" : "PASS"; // Defualt
				format.Set("GT", "1/1");
			}
			vcfline.info_.Add("HRun", "HRun=" + 
							  itoa(allvariant[it->second][i].homoRun));

			format.Add("HR", itoa(allvariant[it->second][i].homoRun));
			format.Add("CI", itoa(allvariant[it->second][i].cipos.first) +","+
						   itoa(allvariant[it->second][i].cipos.second));
			format.Add("CE", itoa(allvariant[it->second][i].ciend.first) +","+
						   itoa(allvariant[it->second][i].ciend.second));
			format.Add("TR", allvariant[it->second][i].target.id + "-" +
						   itoa(allvariant[it->second][i].target.start)+ "-" +
						   itoa(allvariant[it->second][i].target.end));
			format.Add("QR", allvariant[it->second][i].query.id + "-" +
						   itoa(allvariant[it->second][i].query.start) + "-" +
						   itoa(allvariant[it->second][i].query.end));
			// Align End position
			format.Add("AE", itoa(allvariant[it->second][i].target.end));
			format.Add("VT", allvariant[it->second][i].type);
			int qvs = allvariant[it->second][i].query.end - 
					  allvariant[it->second][i].query.start; // Do not +1!!
			int tvs = allvariant[it->second][i].target.end - 
					  allvariant[it->second][i].target.start; // Do not +1!!
			int vs = (qvs > 0) ? qvs : tvs; // Should be tvs if is DEL!
			format.Add("VS", itoa(vs));

			string age (".");
			if (allvariant[it->second][i].isSuccessAlign) {
				age  = (allvariant[it->second][i].isGoodReAlign) ? "T," : "F,"; 
				age += char2str(allvariant[it->second][i].strand)     + "," +
					itoa(allvariant[it->second][i].identity[0].first) + "," +
					itoa(allvariant[it->second][i].identity[0].second)+ "," +
					itoa(allvariant[it->second][i].identity[1].first) + "," +
					itoa(allvariant[it->second][i].identity[1].second)+ "," +
					itoa(allvariant[it->second][i].identity[2].first) + "," +
					itoa(allvariant[it->second][i].identity[2].second);
			}
			format.Add("AGE", age);

			int qn = qryfa.Nlength(allvariant[it->second][i].query.id,
						allvariant[it->second][i].query.start > 100 ? 
						allvariant[it->second][i].query.start - 100 : 0, 
						allvariant[it->second][i].query.end + 100);
			double nr = double (qn) / ( vs + 200);
			format.Add("NR",ftoa(nr));
			vcfline.sample_.push_back(format);

			O << vcfline << "\n";
		}
	}
	O.close();
}

///////////////////////////////////////////////////////////////////////////////
unsigned int RegionMin   (vector<Region> &region) {

	if ( region.empty() ) { std::cerr << "[ERROR] Region is empty, when you're calling RegionMin() function.\n"; exit(1); }
	unsigned int pos = region[0].start;
	for ( size_t i(0); i < region.size(); ++i ) { if (region[i].start < pos ) pos = region[i].start; }
	return pos;
}

unsigned int RegionMax   (vector<Region> &region) {

	if ( region.empty() ) { std::cerr << "[ERROR] Region is empty, when you're calling RegionMax() function.\n"; exit(1); }
	unsigned int pos = region[0].end;
	for ( size_t i(0); i < region.size(); ++i ) { if (region[i].end > pos ) pos = region[i].end; }
	return pos;
}

string ReverseAndComplementary ( string &seq ) {

	string tmpstr;
	for ( int i(seq.size() - 1); i >= 0; --i ) {

		switch(seq[i]) {

			case 'a' : tmpstr.push_back( 't' ); break;
			case 'A' : tmpstr.push_back( 'T' ); break;
			case 'c' : tmpstr.push_back( 'g' ); break;
			case 'C' : tmpstr.push_back( 'G' ); break;
			case 'g' : tmpstr.push_back( 'c' ); break;
			case 'G' : tmpstr.push_back( 'C' ); break;
			case 't' : tmpstr.push_back( 'a' ); break;
			case 'T' : tmpstr.push_back( 'A' ); break;
			default  :
				// Not 'A' 'T' 'C' 'G', then do not complementary
				tmpstr.push_back(seq[i]);
		}
	}

	return tmpstr;
}

vector<VarUnit> MergeVarUnit(vector<VarUnit> &VarUnitVector) {
// CAUTION : Merge vector<VarUnit>, but all new merge result just 
// using the same strand of first element: VarUnitVector[0]!!

	int distDelta = 1;
	bool flag(false);

	VarUnit varunit;
	vector<VarUnit> newVector;
	map<string, unsigned long int> tarPrePos, qryPrePos;
	map<string, string> id2seq;
	for ( size_t i(0); i < VarUnitVector.size(); ++i ) {

		string tarId = VarUnitVector[i].target.id;
		string qryId = VarUnitVector[i].query.id ;
		string id  = VarUnitVector[i].target.id + ":" + VarUnitVector[i].query.id;
		string seq = VarUnitVector[i].tarSeq + "-" + VarUnitVector[i].qrySeq;

		if (!tarPrePos.count(tarId) || !qryPrePos.count(qryId) || !id2seq.count(id)) {

			// The light is on => Get the region!
			if (flag) newVector.push_back(varunit);
			varunit    = VarUnitVector[i];
			id2seq[id] = seq;
			flag       = true; // first time
		} else {

			if (tarPrePos[tarId] > VarUnitVector[i].target.start) {
				std::cerr << "[ERROR]Your target hasn't been sorted.\n"; 
				VarUnitVector[i].target.OutErrReg();
				exit(1);
			}
			if (qryPrePos[qryId] > VarUnitVector[i].query.start) {
				std::cerr << "[ERROR]Your query hasn't been  sorted.\n";
				VarUnitVector[i].query.OutErrReg();
				exit(1);
			}

			if (varunit.target.end + distDelta >= VarUnitVector[i].target.start
				&& varunit.query.end + distDelta >= VarUnitVector[i].query.start
				&& id2seq[id] == seq) {
			
				if (VarUnitVector[i].target.end > varunit.target.end)
					varunit.target.end = VarUnitVector[i].target.end;

				if (VarUnitVector[i].query.end > varunit.query.end)
					varunit.query.end = VarUnitVector[i].query.end; 
			} else {

				newVector.push_back(varunit); 
				varunit = VarUnitVector[i];
			}
		}
		tarPrePos[tarId] = VarUnitVector[i].target.start;
		qryPrePos[qryId] = VarUnitVector[i].query.start;
		id2seq[id]       = seq;
	}
	if (flag) newVector.push_back(varunit);

	return newVector;
}

