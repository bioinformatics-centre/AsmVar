/**********************************
 * Author : Shujia Huang
 * Date   : 2014-08-05 15:20:04
 **********************************/

#ifndef __AGEOPTION_H__
#define __AGEOPTION_H__

class AgeOption {

public :

    int gapOpen;
    int gapExtend;
    int match;
    int mismatch;

	int extendVarFlankSzie;

    bool indel, inv, invl, invr, tdup;
    bool both, revcomTarget, revcomQuery;

public :
//--indel --both --match=1 --mismatch=-1 --go=-10 --ge=-1 
    AgeOption() : match(1), mismatch(-1), gapOpen(-10), gapExtend(-1),
				  extendVarFlankSzie(500),
                  both(true) , revcomTarget(false), revcomQuery(false),
                  indel(true), inv(false), invl(false), invr(false), tdup(false)
	{}

};

#endif

