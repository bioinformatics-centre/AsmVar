/**
 *
 * Author : Shujia Huang
 * Date   : 2014-08-18
 *
 **/
#include "vcf.h"

// For Class Header

void VcfHeader::Add(string mark, string id, string num,
					string type, string description) {


	if (mark == "CHROM") return; // Do not add 'CHROM' here!
	if (data_.count(id)) {
        cerr << "[ERROR] There is already exists " << id << "in header\n";
        exit(1);
    }

	string key = "##%s=<ID=%s>" % (mark, id);
	string val = "##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">" % 
				 (mark, id, num, type, description);
    data_[key] = val;

	return;
}

void VcfHeader::Add(string id, string info) {

	if (data_.count(id)) {
		cerr << "[ERROR] There is already exists " << id << "in header\n";
		exit(1);
	}
	data_[id] = info;
	return;
}


// For Class VCF

void VcfInfo::Add(string id, string info) {

	if (data_.count(id)) {
        cerr << "[ERROR] There is already exists " << id << "in INFO\n";
        exit(1);
    }
    data_[id] = info;
    return;
}

