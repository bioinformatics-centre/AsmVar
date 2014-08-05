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

    bool indel, inv, invl, invr, tdup;
    bool both, revcom1, revcom2;

public :
    AgeOption() : match(1), mismatch(-2), gapOpen(-2), gapExtend(-1),
                  both(false) , revcom1(false), revcom2(false),
                  indel(false), inv(false), invl(false), invr(false),
                  tdup(false)
	{}

};

#endif

