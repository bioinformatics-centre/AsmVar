/* Author : Shujia Huang & Siyang Liu
 * Date   : 2013-12-30 10:12:53	Add gzstream to read file
 *          2013-10-18 15:25:18
 *
 * The main part of variant detection program
 * 1-base system to output
 *
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>
#include "utility.h"
#include "Variant.h"
#include "gzstream.h"

using namespace std;

void Usage        ( const char* prog );
void ReadFileList ( const char* filelist,    vector< string   > & infile );

int main ( int argc, char* argv[] ) {

	char c;
	vector< string > infile, tarRef, qryRef;
	string filelist, outFilePrefix;
	string sampleId, referenceId("ALL");
	while ((c = getopt(argc, argv, "i:l:o:t:r:q:s:h")) != -1) {
		switch ( c ){
            case 'i' : infile.push_back(optarg); break; // .maf file
            case 't' : tarRef.push_back(optarg); break; // Target fa in .maf files. [fa format]
            case 'q' : qryRef.push_back(optarg); break; // Query fa in .maf  files. [fa format]
			case 'o' : outFilePrefix = optarg;   break;
			case 'r' : referenceId   = optarg;   break;
            case 'l' : filelist      = optarg;   break; // .maf file list
			case 's' : sampleId      = optarg;   break; // Sample ID
            case 'h' : Usage( argv[0] );
            default  :
                cerr << "\n[ERROR]Unknow option: -" << c << "\n" << endl;
                Usage( argv[0] );
        }
	}
	if (argc == 1 || infile.empty() || outFilePrefix.empty() || qryRef.empty() || tarRef.empty() ) Usage( argv[0]);
	if (sampleId.empty()) { cerr << "[ERROR] You miss the sample ID by the parameter '-s [sampleId]'\n"; exit(1); }
	cerr << "[INFO] Parameter : " << join(argv, argc) << "\n\n";
	if (toupper(referenceId) == "ALL") {
		cerr << "[WARNING]!!! You're going to use all of the chromosome "
			 << "instead of using '-r' to pick a specific one when in the "
			 << "AGE_Realigne process. That will cause a long time!!! For "
			 << "human genome, it may cause you 7-10 days!\n\n";
	}

	if (!filelist.empty()) ReadFileList(filelist.c_str(), infile);

	Variant variant;
	variant.AssignSample(sampleId); // Assigne the name of sample into 'variant'
	for (size_t i(0); i < qryRef.size(); ++i) { variant.qryfa.Load(qryRef[i]); } // Load query sequence
	for (size_t i(0); i < tarRef.size(); ++i) { variant.tarfa.Load(tarRef[i]); } // Load query sequence

	for (size_t i (0); i < infile.size(); ++i) {

		igzstream I(infile[i].c_str());
        if (!I) {
            cerr << "Cannot open file : " << i + 1 << "\t" << infile[i] << endl;
            exit(1);
        }
        cerr << "***#    Reading file : " << i + 1 << "\t" << infile[i] << "    #*** " << endl;
		string tmp;
        vector<string> vect; 
		while (1) {
		/*  
            # Header ...
            # batch 1
            a score=221 mismap=1e-10
            s chr1        17032965 221 + 35477943 AGTTCTAAGGGCTCCAGTGTACACACATTGCAGAAAC
            s scaffold665   120623 221 -   130124 AGTTCTAAGGGCTCCAGTGTACACACATTGCAGAAAC
            p                                     }}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
            p                                     2567777777777888888889999999:::::;;;;

		*/
			I >> tmp;
			if (I.eof()) {
                break;
            } else if (I.fail()) {
                cerr << "File read ERROR [I/O]!!\n" << tmp << endl;
                exit(1);
            }
			if (tmp[0] == '#') { getline (I, tmp, '\n'); continue; }
			if (tmp[0] != 'a') { 
				cerr << "[ERROR] Format crash. It's not the begin of a MAF alignment block\t" << tmp << endl; 
				exit(1);
			}
            //a score=221 mismap=1e-10
            I >> tmp; split("=", tmp, vect); variant.score  = atoi(vect[1].c_str());  // score=221
            I >> tmp; split("=", tmp, vect); variant.mismap = atof(vect[1].c_str());  // mismap=1e-10
            getline (I, tmp, '\n');
            I >> tmp >> variant.target.id >> variant.target.start >> variant.target.end //Here variant.target.end is just the mapping length
              >> tmp >> tmp >> variant.tarSeq; 
            getline (I, tmp, '\n');
// Just for urgency situation
if (variant.target.id == "M") variant.target.id = "MT"; // Should be deleted after this time

            variant.target.end += variant.target.start; // Now variant.target.end is the end of mapping region
            ++variant.target.start;                     // variant.target.start is 0-base , shift to 1-base

            I >> tmp >> variant.query.id  >> variant.query.start  >> variant.query.end
              >> variant.strand >> tmp    >> variant.qrySeq; getline (I, tmp, '\n');
            variant.query.end += variant.query.start;
            ++variant.query.start;

            do {
				getline (I, tmp, '\n');
			} while (tmp[0] == 'p'); 

			variant.CheckMAF();
			variant.qryfa.CheckFaId(variant.query.id);
			variant.tarfa.CheckFaId(variant.target.id);
			variant.CallHomoRef(referenceId);
			variant.CallnSeq(referenceId);
			variant.CallSNP (referenceId);
			variant.CallInsertion(referenceId); // Don't covert the coordinate which map to the '-' strand here. I'll covert it when calling indel.
			variant.CallDeletion (referenceId); // Don't covert the coordinate which map to the '-' strand here. I'll covert it when calling indel.
			variant.GetMapReg    (); // The coordinate coversion events will be happen in this memerber function!
		}
		I.close();
	}
	variant.CallSV(); // It's SV not Indel, which cannot be called by a single alignment. The most important part of these program!!
	variant.CallClipReg();
	variant.CallNomadic();
	variant.Filter();      // Filter the indels' regions which in nosolution regions. Maybe we don't need it! Maybe I should just call this function in AGE_Realign()
	cerr << "[INFO] Doing re-aligne process ...\n";
	variant.AGE_Realign(referenceId);

	cerr << "[INFO] Outputting information into files.\n";
	variant.Output   (outFilePrefix + ".svd"    );
	cerr << "      -- " << outFilePrefix + ".svd\n";
	variant.OutputGap(outFilePrefix + ".gap.bed");
	cerr << "      -- " << outFilePrefix + ".gap.bed\n";
	variant.Output2VCF(referenceId, outFilePrefix + ".vcf");
	cerr << "      -- " << outFilePrefix + ".vcf\n";
	variant.Summary  (outFilePrefix + ".summary"); 
	cerr << "      -- " << outFilePrefix + ".summary\n";

	cerr << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> All Done <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

	return 0;
}

void ReadFileList(const char* filelist, vector< string > & infile) {

    igzstream I(filelist);
    if (!I) { cerr << "Cannot open file : " << filelist << endl; exit(1); }
    string tmp;
    while (1) {

        I >> tmp;
        if (I.eof()) break;

        infile.push_back(tmp);
        getline(I, tmp, '\n');
    }
    I.close();

    return;
}

void Usage(const char* prog) {

    cerr << "Version : 0.0.0 ( 2013-10-18 )                                                          \n"
        << "Author  : Shujia Huang                                                                   \n"
        << "Created : 2013-10-18                                                                     \n"
        << "Last Modify : 2013-11-19 17:31:49 Fix bugs                                               \n"
		<< "\nCaution: The the input axt file should be the regular axt format. It means do not covert\n"
		<< "           the coordinates of querys which align to '-' strand youself. It's not right and excrescent!\n"
        << "\nUsage : " << prog << " [Options] ""-i [axt infile] -o output                           \n"
        << "\n   Options  :                                                                          \n"
        << "       -i  [str]   maf file. require!                                                    \n"
        << "       -l  [str]   maf file list. [none]                                                 \n"
        << "       -o  [str]   Output file prefix. require!                                          \n"
        << "       -t  [str]   target sequence, fa format. require!                                  \n"
        << "       -q  [str]   Query  sequence, fa format. require!                                  \n"
        << "       -r  [str]   Specific chromosome ID. [ALL]                                         \n"
		<< "       -s  [str]   Sample ID. require!                                                   \n"
        << "       -h          Output this help information.                                         \n"
        << endl;

    exit(1);
}

