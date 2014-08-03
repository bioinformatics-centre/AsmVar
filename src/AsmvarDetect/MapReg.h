/**************************
 * Author : Shujia Huang
 * Date   : 2014-08-01
 **************************/
#ifndef __MAPREG_H__
#define __MAPREG_H__

#include <iostream>
#include "Region.h"

class MapReg { // Mapping regions

public:
    Region target; // target or said reference
    Region query;  // the mapping one
    char  strand;  // Mapping strand

    long score;
    double mismap; // mismatch probability

public:
    void OutErrReg() {
        std::cerr << "# " << score << "\t" << mismap << "\t" << target.id 
				  << "\t" << target.start  << "\t"   << target.end << "\t"
                  << query.id << "\t"  << query.start << "\t" << query.end  
				  << "\t"    << strand << endl;
    }
};

#endif
