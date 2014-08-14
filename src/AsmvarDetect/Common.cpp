/*
 * Author : Shujia Huang
 * Date   : 2014-08-15
 * 
 */
#include "Common.h"

unsigned int RegionMin   ( vector<Region> & region ) {

    if ( region.empty() ) { std::cerr << "[ERROR] Region is empty, when you're calling RegionMin() function.\n"; exit(1); }
    unsigned int pos = region[0].start;
    for ( size_t i(0); i < region.size(); ++i ) { if (region[i].start < pos ) pos = region[i].start; }
    return pos;
}

unsigned int RegionMax   ( vector<Region> & region ) {

    if ( region.empty() ) { std::cerr << "[ERROR] Region is empty, when you're calling RegionMax() function.\n"; exit(1); }
    unsigned int pos = region[0].end;
    for ( size_t i(0); i < region.size(); ++i ) { if (region[i].end > pos ) pos = region[i].end; }
    return pos;
}

string ReverseAndComplementary ( string & seq ) {

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

vector<VarUnit> MergeVarUnit( vector<VarUnit>& VarUnitVector ){
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


