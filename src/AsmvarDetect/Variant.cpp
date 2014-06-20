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

// Used to calculate the 'n' length in class 'VarUnit' and 'Variant'
unsigned int NLength ( string & str ) {

	uint len(0);
	for ( size_t i(0); i < str.size(); ++i )
		if ( str[i] == 'N' || str[i] == 'n' ) ++len;

	return len;
}

void Variant::CallSNP () {

	assert ( tarSeq.length() == qrySeq.length() );

	VarUnit tmpgap;
    tmpgap.target.id = target.id;
    tmpgap.query.id  = query.id;
    tmpgap.strand    = strand;
	tmpgap.score     = score;   // 2014-06-20 09:58:32
	tmpgap.mismap    = mismap;  // 2014-06-20 09:58:39
    tmpgap.type      = "SNP";

	int tarPos(target.start), qryPos(query.start);
	int tarI(0), qryI(0);
	for ( int i(0); i < tarSeq.length(); ++i ) {
	// e.g. : tarSeq = "-ab-c-t--",
	//        qrySeq = "aa-ccdcdt",
		tarPos += tarI;
		qryPos += qryI;
		tarI = tarSeq[i] == '-' ? 0 : 1;
		qryI = qrySeq[i] == '-' ? 0 : 1;
		if (tarI == 0 || qryI  == 0) continue;
		if (toupper(tarSeq[i]) == 'N' || toupper(qrySeq[i]) == 'N') continue;

		if ( toupper(tarSeq[i]) != toupper(qrySeq[i]) ) {

			tmpgap.target.start = tarPos; tmpgap.target.end = tarPos;
			tmpgap.query.start  = qryPos; tmpgap.query.end  = qryPos;
			tmpgap.ConvQryCoordinate( qryfa.fa[query.id].length() ); // coordinates uniform to the positive strand!
			snp.push_back( tmpgap );
		}
	}
}

void Variant::CallInsertion () { 
// Actually, we just have to call  the target gap regions.	
// All the coordinate of query should be uniform to the positive strand!
	vector< VarUnit > gap = CallGap ( target, tarSeq, query, qrySeq, strand, score, mismap, "Ins" );
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
	vector< VarUnit > gap = CallGap ( query, qrySeq, target, tarSeq, strand, score, mismap, "Del" );
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

	// Call the query coverting function here to make the '-' strand coordinates of query be the same as '+' strand!
	ConvQryCoordinate(); // ConvQryCoordinate() is a memerber function of class 'MAF'

	MapReg mpI;
	mpI.target = target; mpI.query = query; mpI.strand = strand;
	mpI.score  = score ; mpI.mismap= mismap; // 2014-06-20 10:39:54
	mapreg[query.id].push_back( mpI );       // The 'key' is query.id

	target.info = query.id;
	query.info  = target.id;
	maptar[target.id].push_back( target ); // stored the mapped target regions here
	mapqry[query.id].push_back ( query  ); // stored the mapped query  regions here
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

void Variant::CallSV () { // Call Stuctural variants, not indels!!! 
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
				} else if( (labs(qry2tarOrder[pos2] - qry2tarOrder[pos1]) == 1) && (labs(qry2tarOrder[pos] - qry2tarOrder[pos2]) == 1) ){ 

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
				cerr << "\n[ERROR] Program bugs! " << endl;
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

	if ( qryfa.fa.empty() ) { cerr << "[WARNING] No Clip regions. Because the query fa is empty!" << endl; return; }
	VarUnit tmp;
	unsigned int rStart, rEnd;
	for ( map< string, vector<Region> >::iterator it( mapqry.begin() ); it != mapqry.end(); ++it ) {

		rStart = RegionMin( it->second );
		rEnd   = RegionMax( it->second );

		if ( rStart == 0 || rEnd == 0 ) { 

			cerr << "rStart == 0 || rEnd == 0" << "\nrStart: " << rStart << "\trEnd: " << rEnd << endl;
			for ( size_t i(0); i < it->second.size(); ++i )	
				cerr << "Size: " << it->second.size() << "\tit->second: " << it->second[i].id 
					 << "\t" << it->second[i].start   << "\t" << it->second[i].end << endl;
			exit(1);
		}

		if ( rStart == 1 && rEnd == qryfa.fa[it->first].length() ) summary["Full-Align"]++;

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

	if ( qryfa.fa.empty() ) { cerr << "[WARNING] No Nomadic regions. Because the query fa is empty!" << endl; return; }
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

	map< string,vector<Region> > mapNosolution;
	map<string, size_t> index;

	for ( size_t i(0); i < nosolution.size(); ++i ) mapNosolution[nosolution[i].query.id].push_back( nosolution[i].query );
	for ( map< string,vector<Region> >::iterator it( mapNosolution.begin() ); it != mapNosolution.end(); ++it ) {
		sort ( mapNosolution[it->first].begin(), mapNosolution[it->first].end(), SortRegion );
		index[it->first] = 0;
	}

	FilterReg ( mapNosolution, index, insertion ); // Filter insertion regions which in mapNosolution
	FilterReg ( mapNosolution, index, deletion  ); // Filter deletion  regions which in mapNosolution
	FilterReg ( mapNosolution, index, snp       ); // Filter SNP which in mapNosolution

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

void Variant::FilterReg( map< string,vector<Region> > tarregion, map<string, size_t> index, vector<VarUnit>& region ) {
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
    if ( !O ) { cerr << "Cannot write to file : " << file << endl; exit(1); }
	for ( map<string, unsigned int>::iterator pt( summary.begin() ); pt!= summary.end(); ++pt ) 
		O << pt->first << "\t" << pt->second << "\n";
	
	O << "qryCovlength/querylength  " << double ( summary["qryCovlength"] ) / qryfa.length << "\n";
	O << "tarCovlength/targetlength " << double ( summary["tarCovlength"] ) / tarfa.length << "\n";
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
    if ( !O ) { cerr << "Cannot write to file : " << file << endl; exit(1); }

	Output ( insertion, O ); //
	Output ( deletion,  O ); //
	Output ( inversion, O );
	Output ( simulreg,  O ); //
	Output ( clipreg,   O ); //
	Output ( nomadic,       O ); //
	Output ( nosolution,    O ); 
	Output ( translocation, O );

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
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n";
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
		if ( !isNext ) {
			O   << snp[i].target.id << "\t" << snp[i].target.start << "\t.\t" << snp[i].tarSeq << "\t" << snp[i].qrySeq << "\t255\tPASS\t.\tGT:ST:VT:QR:MS:MIP\t"
				<< "./.:" << snp[i].strand  << ":" + snp[i].type + ":" + snp[i].query.id + "-" + itoa(snp[i].query.start) + "-" + itoa(snp[i].query.end) + ":" 
                   + itoa( snp[i].score ) + ":" << snp[i].mismap << "\n";
		} else {
			cerr << "#\t"  << snp[i].target.id << "\t" << snp[i].target.start << "\t.\t" << snp[i].tarSeq << "\t" << snp[i].qrySeq << "\t255\tPASS\t.\tGT:ST:VT:QR:MS:MIP\t"
                 << "./.:" <<  snp[i].strand   <<  ":" + snp[i].type + ":" + snp[i].query.id + "-" + itoa(snp[i].query.start) + "-" + itoa(snp[i].query.end) + ":" 
                   + itoa( snp[i].score ) + ":" << snp[i].mismap << "\n";
		}
    }
	O.close();
}

void Variant::OutputGap( string file ) {
// The inter-gaps between different scaffolds in the same chromosome need to calculate.
// Abort : 2013-11-04 13:20:30 !!!

	ofstream O ( file.c_str() );
	if ( !O ) { cerr << "Cannot write to file : " << file << endl; exit(1); }

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
                                     string type )    // insertion or deletion
{
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
		gap.target.start = left.target.end;   // Should +1? what if the size larger than the size of query after +1. so I don't want to +1!
		gap.target.end   = right.target.start;// I don't want to -1 ! If left.query overlap with right.query, we will find gap.query.start >= gap.query.end
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

unsigned int RegionMin   ( vector<Region> & region ) {

	if ( region.empty() ) { cerr << "[ERROR] Region is empty, when you're calling RegionMin() function.\n"; exit(1); }
	unsigned int pos = region[0].start;
	for ( size_t i(0); i < region.size(); ++i ) { if (region[i].start < pos ) pos = region[i].start; }
	return pos;
}
unsigned int RegionMax   ( vector<Region> & region ) {

	if ( region.empty() ) { cerr << "[ERROR] Region is empty, when you're calling RegionMax() function.\n"; exit(1); }
	unsigned int pos = region[0].end;
	for ( size_t i(0); i < region.size(); ++i ) { if (region[i].end > pos ) pos = region[i].end; }
	return pos;
}

string ReverseAndComplementary ( string & seq ) {

	string tmpstr;
	for ( int i(seq.size() - 1); i >= 0; --i ) {

		if ( toupper( seq[i] ) == 'A' ) {

			if ( seq[i] == 'a' ) tmpstr.push_back( 't' );
			if ( seq[i] == 'A' ) tmpstr.push_back( 'T' );
		} else if ( toupper( seq[i] ) == 'C' ) {
			if ( seq[i] == 'c' ) tmpstr.push_back( 'g' );
			if ( seq[i] == 'C' ) tmpstr.push_back( 'G' );
		} else if ( toupper( seq[i] ) == 'G' ) {
			if ( seq[i] == 'g' ) tmpstr.push_back( 'c' );
			if ( seq[i] == 'G' ) tmpstr.push_back( 'C' );
		} else if ( toupper( seq[i] ) == 'T' ) {
			if ( seq[i] == 't' ) tmpstr.push_back( 'a' );
			if ( seq[i] == 'T' ) tmpstr.push_back( 'A' );
		} else {
			tmpstr.push_back( seq[i] );
		}
	}

	return tmpstr;
}




