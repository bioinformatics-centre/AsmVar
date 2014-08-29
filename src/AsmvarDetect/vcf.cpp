/**
 *
 * Author : Shujia Huang
 * Date   : 2014-08-18
 *
 **/
#include "vcf.h"

// For Class Header
void VcfHeader::DefualtHeader() {

	Add("FORMAT", "GT", "1", "String", "Genotype");
	Add("FORMAT", "AE", "1", "Integer", "End position of reference call block");
	Add("FORMAT", "AGE", "1", "String", 
		"AGE aligment information. (T/F,Strand,ave_base,ave_iden,left_base,left_iden,right_base,right_iden)");
	Add("FORMAT", "CI", "2", "Integer", 
		"Confidence interval around POS for imprecise variants");
	Add("FORMAT", "CE", "2", "Integer", 
        "Confidence interval around END for imprecise variants");
	Add("FORMAT", "HR", "1", "Integer", "Homo-Run of reference");
	Add("FORMAT", "MS", "1", "Float", "Mismatch probability In LAST Align");
	Add("FORMAT", "NR", "1", "Float", 
		"The 'N' ratio around variant region on query");
	Add("FORMAT", "SC", "1", "Integer", "Alignment score In LAST Align");
	Add("FORMAT", "TR", "1", "String",  "Variant regions on reference");
	Add("FORMAT", "QR", "1", "String",  "Variant regions on query");
	Add("FORMAT", "VS", "1", "Integer", "SV Size");
	Add("FORMAT", "VT", "1", "String" , "SV Type");

	Add("INFO", "HRun", "1", "Integer", "Homo-Run of reference");
	Add("##FILTER=<ID=INTERGAP>", 
		"##FILTER=<ID=INTERGAP,Description=\"Unmapped reference region.\">");
	Add("##FILTER=<ID=NCALL>", 
		"##FILTER=<ID=NCALL,Description=\"Covered by Query, but it's 'N' on Reference or Query or both of them are 'N'\">");
	Add("##FILTER=<ID=REFCALL>", 
		"##FILTER=<ID=REFCALL,Description=\"Represents a homozygous reference call\">");
	Add("##FILTER=<ID=AGEFALSE>","##FILTER=<ID=AGEFALSE,Description=\"Filtered in AGE realignment process\">");

	return;
}

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
        cerr << "[ERROR] There is already exists " << id << "in FORMAT. ";
		cerr << "You can use Set(string, string) to reset the data.\n";
        exit(1);
    }
    data_[id] = dat;
    return;
}

void VcfFormat::Set(string id, string dat) {

	if (!data_.count(id)) {
        cerr << "[ERROR] " << id << " is still not exists in FORMAT. ";
		cerr << "You should call Add() first.\n";
        exit(1);
    }
    data_[id] = dat;
    return;
}

string VcfFormat::GetFormat() {

	string format = "GT";
	for (map<string, string>::iterator it(data_.begin());
		 it != data_.end();
         ++it) {
		if (it->first != "GT") format += ":" + it->first;
	}
	return format;
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
		   info_.Combine();
	if (sample_.size() > 0) { 
		format_ = sample_[0].GetFormat();
		data   += "\t" + format_;
	}
	for (size_t i(0); i < sample_.size(); ++i) 
		data += "\t" + sample_[i].Combine();

	return data;
}

ostream &operator << (ostream &o, VCF vcf) {

	o << vcf.Combine();
    return o;
}




