/**********************************
 * Author : Shujia Huang
 * Date   : 2014-08-05 22:49:56
 *********************************/

#ifndef __AGEALIGNMENT_H__
#define __AGEALIGNMENT_H__

//--- AGE includes ---
#include "AGEaligner.h"
#include "Sequence.h"
#include "Scorer.h"
#include "AgeOption.h"

#ifdef AGE_TIME
#include <sys/time.h>
#endif

//--- Other include ---
#include "VarUnit.h"

class AgeAlignment {

public:

	AgeAlignment();
	bool Init(VarUnit &v, AgeOption op);
	bool Align(string &tarFa, string &qryFa);

private:

	void ExtendVariant(unsigned long int, unsigned long int, int extandFlankSzie);
	bool IsHugeMemory(unsigned long int n, unsigned long int m);

private:

	VarUnit   vu;
	AgeOption opt;
	bool isInit;
};

#endif


