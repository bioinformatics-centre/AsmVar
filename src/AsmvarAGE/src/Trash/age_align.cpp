/****************************************************************************
 *
 *     AGE -- Alignment with Gap Excision
 *     Copyright (C)  Alexej Abyzov
 *                
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Creative Commons license
 * (Attribution-NonCommerical).
 * See terms at http://creativecommons.org/licenses/by-nc/2.5/legalcode
 *
 * Author: Alexej Abyzov
 */

//--- C/C++ includes ---
#include <cstring>
#include <cstdlib>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#ifdef AGE_TIME
#include <sys/time.h>
#endif
using namespace std;

//--- AGE includes ---
#include "AGEaligner.hh"
#include "Sequence.hh"
#include "Scorer.hh"

bool parseValue(string &str,short &val)
{
	if (str == "-") return false;

	int n = str.length(),i = 0;
	string s("");
	if (str[i] == '-') s += str[i++];
	while (i < n && isdigit(str[i])) s += str[i++];
	if (i != n)   return false;
	val = strtol(s.c_str(),NULL,10);
	return true;
}

bool parseStartEnd(string &str,int &start,int &end)
{
	if (str == "-") return false;

	int n = str.length(),i = 0;
	string s(""),e("");
	while (i < n && isdigit(str[i])) s += str[i++];
	if (i == n)          return false;
	if (str[i++] != '-') return false;
	while (i < n && isdigit(str[i])) e += str[i++];
	if (i != n)          return false;

	int sv = -1,ev = -1;
	if (s.length() > 0) sv = strtol(s.c_str(),NULL,10);
	if (e.length() > 0) ev = strtol(e.c_str(),NULL,10);
	if (sv == 0 || ev == 0)          return false;
	if (sv > 0 && ev > 0 && ev < sv) return false;
	start = sv;
	end   = ev;
	return true;
}

int parseArgs(string line,string *args,int max_args)
{
	int ret = 0, i = 0;
	while (i < line.length()) {
		args[ret] = "";
		while (i < line.length() && line.at(i) == ' ') i++; // Skipping spaces
		while (i < line.length() && line.at(i) != ' ')
			args[ret] += line.at(i++);
		if (args[ret].length() > 0) ret++;
	}

	return ret;
}

int main(int argc,char *argv[])
{
	const static short DEFAULT_MATCH      =  1;
	const static short DEFAULT_MISMATCH   = -2;
	const static short DEFAULT_GAP_OPEN   = -2;
	const static short DEFAULT_GAP_EXTEND = -1;

	int flag = 0;
	stringstream sout;

	string exe = argv[0];
	exe = exe.substr(exe.rfind('/') + 1);
	string usage = "Usage:\n\t";
	usage += exe; usage += "\n";
	sout.str("");sout<<DEFAULT_MATCH;
	usage += "\t\t[-version]\n";
	usage += "\t\t[-indel|-tdup|-inv|-invl|-invr]\n";
	usage += "\t\t[-match=value:"; usage += sout.str(); usage += "]";
	sout.str("");sout<<DEFAULT_MISMATCH;
	usage += " [-mismatch=value:"; usage += sout.str(); usage += "]\n";
	sout.str("");sout<<DEFAULT_GAP_OPEN;
	usage += "\t\t[-go=value:";    usage += sout.str(); usage += "]";
	sout.str("");sout<<DEFAULT_GAP_EXTEND;
	usage += " [-ge=value:";       usage += sout.str(); usage += "]\n";
	usage += "\t\t[-both] [-revcom1] [-revcom2]\n";
	usage += "\t\t[-coor1=start-end] [-coor2=start-end]";
	usage += " file1.fa file2.fa";
	string age_version = "AGE ";
#ifdef AGE_VERSION
	age_version += AGE_VERSION;
#else
	age_version += "v???";
#endif
	age_version += "\n";

	cin.seekg(0,ios::end);
	int cin_length = cin.tellg();
	cin.seekg(0,ios::beg);

	if (argc == 2) {
		string tmp(argv[1]);
		if (tmp == "-version") {
			cout<<age_version<<endl;
			return 0;
		}
	}

	if (argc < 3 && cin_length <= 0) {
		cerr<<usage<<endl;
		return 0;
	}

	const static int MAX_ARGS = 1024;
	string args[MAX_ARGS],line;
	int    n_args = 0;

	if (argc > 1) // If we have input arguments
		for (int i = 1;i < argc;i++) {
			args[n_args] = argv[i];
			line += " ";
			line += args[n_args++];
			if (n_args == MAX_ARGS) {
				cerr<<"Too many input arguments. Skipping!"<<endl;
				break;
			}
		}
	else while (cin_length > 0) { // Taking lines from redirection
		getline(cin,line);
		cin_length -= line.length() + 1;
		if (line.length() > 0) {
			n_args = parseArgs(line,args,MAX_ARGS);
			break;
		}
	}

	string    old_file1(""),   old_file2("");
	Sequence *old_seq1 = NULL,*old_seq2 = NULL;
	while (n_args > 0) {
		short match      = DEFAULT_MATCH;
		short mismatch   = DEFAULT_MISMATCH;
		short gap_open   = DEFAULT_GAP_OPEN;
		short gap_extend = DEFAULT_GAP_EXTEND;
		string file1(""),file2("");
		int start1 = -1,end1 = -1;
		int start2 = -1,end2 = -1;
		bool revcom1 = false, revcom2 = false, both = false, error = false;
		bool version = false;
		for (int i = 0;i < n_args;i++) {
			if (args[i].at(0) == '-') { // Options
				string option = args[i];
				string begin = option.substr(0,7);
				if (begin == "-match=") {
					string value = option.substr(7);
					if (!parseValue(value,match)) {
						cerr<<line<<endl;
						cerr<<"Error in match specification '"<<option<<"'."<<endl;
						cerr<<"Using default value "<<match<<"."<<endl;
					}
				} else if (option.substr(0,10) == "-mismatch=") {
					string value = option.substr(10);
					if (!parseValue(value,mismatch)) {
						cerr<<line<<endl;
						cerr<<"Error in mismatch specification '"<<option<<"'."<<endl;
						cerr<<"Using default value "<<mismatch<<"."<<endl;
					}
				} else if (option.substr(0,4) == "-go=") {
					string value = option.substr(4);
					if (!parseValue(value,gap_open)) {
						cerr<<line<<endl;
						cerr<<"Error in gap open specification '"<<option<<"'."<<endl;
						cerr<<"Using default value "<<gap_open<<"."<<endl;
					}
				} else if (option.substr(0,4) == "-ge=") {
					string value = option.substr(4);
					if (!parseValue(value,gap_extend)) {
						cerr<<line<<endl;
						cerr<<"Error in gap extend specification '"<<option<<"'."<<endl;
						cerr<<"Using default value "<<gap_extend<<"."<<endl;
					}
				} else if (begin == "-coor1=") {
					string end = option.substr(7);
					if (!parseStartEnd(end,start1,end1)) {
						cerr<<line<<endl;
						cerr<<"Error in coordinate specification '"<<option<<"'."<<endl;
						cerr<<usage<<endl;
						error = true;
						break;
					}
				} else if (begin == "-coor2=") {
					string end = option.substr(7);
					if (!parseStartEnd(end,start2,end2)) {
						cerr<<line<<endl;
						cerr<<"Error in coordinate specification '"<<option<<"'."<<endl;
						cerr<<usage<<endl;
						error = true;
						break;
					}
				} else if (option == "-version") {
					version = true;
				} else if (option == "-indel") {
					flag |= AGEaligner::INDEL_FLAG;
				} else if (option == "-inv") {
					flag |= AGEaligner::INVERSION_FLAG;
				} else if (option == "-invl") {
					flag |= AGEaligner::INVL_FLAG;
				} else if (option == "-invr") {
					flag |= AGEaligner::INVR_FLAG;
				} else if (option == "-tdup") {
					flag |= AGEaligner::TDUPLICATION_FLAG;
				} else if (option == "-revcom1") {
					revcom1 = true;
				} else if (option == "-revcom2") {
					revcom2 = true;
				} else if (option == "-both") {
					both = true;
				} else {
					cerr<<line<<endl;
					cerr<<"Unknown option '"<<option<<"'."<<endl;
					cerr<<usage<<endl;
					error = true;
					break;
				}
				continue;
			}

			if (file1.length()      <= 0) file1 = args[i];
			else if (file2.length() <= 0) file2 = args[i];
			else {
				cerr<<line<<endl;
				cerr<<"Too many input files '"<<args[i]<<"'."<<endl;
				cerr<<usage<<endl;
				error = true;
				break;
			}
		}

		if (version) {
			cout<<age_version<<endl;
			return 0;
		}

		if (flag == 0) flag = AGEaligner::INDEL_FLAG;

		int n_modes = 0;
		if (flag & AGEaligner::INDEL_FLAG)        n_modes++;
		if (flag & AGEaligner::TDUPLICATION_FLAG) n_modes++;
		if (flag & AGEaligner::INVR_FLAG)         n_modes++;
		if (flag & AGEaligner::INVL_FLAG)         n_modes++;
		if (flag & AGEaligner::INVERSION_FLAG)    n_modes++;

		if (n_modes != 1) {
			cerr<<"Error in mode specification. ";
			if (n_modes == 0)
				cerr<<"No mode is specified."<<endl;
			if (n_modes > 1)
				cerr<<"More than one mode is specified."<<endl;
			error = true;
		}

		if (file1 != "" && file2 == "") cerr<<"No second file given."<<endl;

		if (!error && file1 != "" && file2 != "") {


			if (file1 != old_file1) {
				Sequence::deleteSequences(old_seq1);
				//old_seq1  = new Sequence(file1);
				old_seq1  = Sequence::parseSequences(file1);
				old_file1 = file1;
			}

			if (file2 != old_file2) {
				Sequence::deleteSequences(old_seq2);
				//old_seq2  = new Sequence(file2);
				old_seq2  = Sequence::parseSequences(file2);
				old_file2 = file2;
			}

			for (Sequence *s1 = old_seq1;s1;s1 = s1->next())
				for (Sequence *s2 = old_seq2;s2;s2 = s2->next()) {

					Sequence *seq1 = s1->substr(start1,end1);
					Sequence *seq2 = s2->substr(start2,end2);

					if (revcom1) seq1->revcom();
					if (revcom2) seq2->revcom();

					// Sequence seq1("ACTGGTGTCAACTG"),seq2("ACTGACTG");
					// Sequence seq1("GGACTGGACTG"),   seq2("CACTGGACTG");
					// Sequence seq1("CCCCAGCTTT"),    seq2("AGCCTTT");
					// Sequence seq1("AGCCTTT"),       seq2("AGCCTTT");

#ifdef AGE_TIME
					timeval ali_s,ali_e;
					gettimeofday(&ali_s,NULL);
#endif

					Scorer scr(match,mismatch,gap_open,gap_extend);
					if (both) {
						Sequence *seq3 = seq2->clone();
						seq3->revcom();
						AGEaligner aligner1(*seq1,*seq2);
						AGEaligner aligner2(*seq1,*seq3);
						bool res1 = aligner1.align(scr,flag);
						bool res2 = aligner2.align(scr,flag);
						if (!res1 && !res2) cerr<<"No alignment made."<<endl;
						else if (aligner1.score() >= aligner2.score())
							aligner1.printAlignment();
						else 
							aligner2.printAlignment();
						delete seq3;
					} else {
						AGEaligner aligner(*seq1,*seq2);
						if (aligner.align(scr,flag)) aligner.printAlignment();
						else cerr<<"No alignment made."<<endl;
					}

#ifdef AGE_TIME
					gettimeofday(&ali_e,NULL);
					cout<<endl<<"Alignment time is "<<ali_e.tv_sec - ali_s.tv_sec +
						(ali_e.tv_usec - ali_s.tv_usec)/1e+6<<" s"<<endl<<endl;
#endif

				}
		}

		n_args = 0;
		while (cin_length > 0) {
			getline(cin,line);
			cin_length -= line.length() + 1;
			if (line.length() > 0) {
				n_args = parseArgs(line,args,MAX_ARGS);
				break;
			}
		}
	}

	Sequence::deleteSequences(old_seq1);
	Sequence::deleteSequences(old_seq2);

	return 1;
}
