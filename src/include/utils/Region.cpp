/**
 *
 * Author : Shujia Huang
 * Date   : 2014-08-21
 *
 **/

#include "Region.h"

vector<Region> MergeRegion(vector<Region> &regVect, int delta=1) {
// CAUTION: The 'info' value in each member of regVect will been clean
// after merge!

	vector<Region> newVect;
	map<string, long int> prePos;

	Region reg;
	bool flag(false);
	for (size_t i(0); i < regVect.size(); ++i) {

		if (regVect[i].start > regVect[i].end) {
			cerr << "[ERROR] MergeRegion ERROR! Your region start > end, "
				 << "which does not allow when calling MergeRegion()!\n";
			exit(1);
		}
		
		if (prePos.count(regVect[i].id)) {

			if (prePos[regVect[i].id] > regVect[i].start) {
				cerr << "[ERROR]Your target hasn't been sorted.\n";
				regVect[i].OutErrReg();
				exit(1);
			}

			if (reg.end + delta >= regVect[i].start) {
				if (reg.end < regVect[i].end) reg.end = regVect[i].end;
			} else {
				newVect.push_back(reg);
				reg = regVect[i];
			}
		} else {

			if (flag) newVect.push_back(reg);
			reg  = regVect[i];
			flag = true;
		}
		prePos[regVect[i].id] = regVect[i].start;
	}
	if (flag) newVect.push_back(reg);

	return newVect;
}

