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

	string key = "##" + mark + "=<ID=" + id + ">";
	string val = "##" + mark + "=<ID=" + id + ",Number="+ num + 
				 ",Type=" + type + ",Description=\""+ description + "\">";
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


// For Class VcfInfo
void VcfInfo::Add(string id, string info) {

	if (data_.count(id)) {
        cerr << "[ERROR] There is already exists " << id << "in INFO\n";
        exit(1);
    }
    data_[id] = info;
    return;
}

string VcfInfo::Combine() { 

	string info;
	for (map<string, string>::iterator it(data_.begin());
         it != data_.end();
         ++it) {

		if (it == data_.begin()) {
			info = it->second;
		} else {
			info += ":" + it->second;
		}
	}

	return info;
}


// For class VcfFormat
void VcfFormat::Add(string id, string dat) {

    if (data_.count(id)) {
        cerr << "[ERROR] There is already exists " << id << "in FORMAT\n";
        exit(1);
    }
    data_[id] = dat;
    return;
}

string VcfFormat::Combine() {

	string gt = data_["GT"];
	data_.erase("GT");
	
	string format = gt;
	for (map<string, string>::iterator it(data_.begin()); 
		 it != data_.end(); 
		 ++it) format += ":" + it->second;

	data_["GT"] = gt; // Set back!
	return format;
}


// VCF class
string VCF::Combine() {

	string data;
	data = chrom_     + "\t" + 
		   itoa(pos_) + "\t" + 
		   Id_        + "\t" +
		   ref_       + "\t" +
           alt_       + "\t" +
           itoa(qual_)+ "\t" +
           filters_   + "\t" +
		   info_.Combine() + "\t" + format_;
	for (size_t i(0); i < sample_.size(); ++i) 
		data += "\t" + sample_[i].Combine();

	return data;
}

ostream &operator << (ostream &o, VCF vcf) {

	o << vcf.Combine();
    return o;
}




