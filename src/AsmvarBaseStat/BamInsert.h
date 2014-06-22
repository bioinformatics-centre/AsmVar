/*
 * Author : Shujia Huang & Siyang Liu
 * Date   : 2013-11-22 14:37:13
 *
 * Caculate insertsize by different @RG
 *
 */

#ifndef BAMINSERT_H
#define BAMINSERT_H

#include <iostream>
#include <fstream>
#include <utility>      // std::pair
#include <map>
#include <math.h>       /* sqrt */

#include "api/BamMultiReader.h"
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

// insert size calculate and return a inner average fragment coverage depth
bool GetInsersizeFromBam ( string bamInfile, map< string, pair<int, int> >& is, uint32_t n=200000 ); // RG => pair<int, int>. 'first' is mean, 'second' is SD
pair<int, int> MeanAndStddev ( vector<int> & v ); // pair.first = mean, pair.second = standard deviation

#endif

