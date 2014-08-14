/*
 * Author : Shujia Huang
 * Date   : 2014-08-15
 * 
 * This is a common.h for AsmvarDetect
 * 
 */
#ifndef __COMMON_H__
#define __COMMON_H__

#include <vector>
#include "Region.h"
#include "VarUnit.h"

unsigned int RegionMin( vector<Region> & region );
unsigned int RegionMax( vector<Region> & region );
string ReverseAndComplementary( string & seq );
vector<VarUnit> MergeVarUnit( vector<VarUnit>& VarUnitVector );

#endif

