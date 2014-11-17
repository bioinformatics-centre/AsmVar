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
	Add("FORMAT", "AGE", "1", "String", 
		"AGE aligment information. (T/F,Strand,ave_base,ave_iden,left_base,left_iden,right_base,right_iden)");
	Add("FORMAT", "HR", "1", "Integer", 
		"Largest Contiguous Homopolymer Run of Variant Allele In Either Direction");
	Add("FORMAT", "MS", "1", "Float", "Mismatch probability In LAST Align");
	Add("FORMAT", "NR", "1", "Float", 
		"The 'N' ratio around variant region in query");
	Add("FORMAT", "AS", "1", "Integer", "Alignment score In LAST Aligner");
	Add("FORMAT", "END", "1", "Integer", "Stop position of the interval");
	Add("FORMAT", "TR", "1", "String",  "Variant regions in reference");
	Add("FORMAT", "QR", "1", "String",  "Variant regions in query");
	Add("FORMAT", "VS", "1", "Integer", "SV Size. For SNP it'll be 1");
	Add("FORMAT", "VT", "1", "String" , 
		"SV Type. Including: REFCALL(Mapped, and it's homozygous reference), REFGAP(Mapped, but contained 'N' in REF-Region), INTRAGAP(Mapped, but contained 'N' in ALT-Region), INTERGAP(Unmapped region in reference), SNP, INS, DEL, MNP, REPLACEMENT, INV, TRANS; For TRANS the format is : VT=>TRANS#TRANSLOCATED_TR#TRANSLOCATED_QR. For the variants that the REF and ALT alleles are the same, the ALT alleles will be replaced by 'N' string");

	Add("##FILTER=<ID=INTERGAP>", 
		"##FILTER=<ID=INTERGAP,Description=\"Unmapped reference region.\">");
	Add("##FILTER=<ID=INTRAGAP>", 
		"##FILTER=<ID=INTRAGAP,Description=\"Mapped, but contained N in Query\">");
	Add("##FILTER=<ID=REFGAP>", 
		"##FILTER=<ID=REFGAP,Description=\"Mapped, but contained N in Reference\">");
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

	if (data_.empty()) return ".";

	string info;
	for (map<string, string>::iterator it(data_.begin());
         it != data_.end();
         ++it) {

		if (it == data_.begin()) {
			info = it->second;
		} else {
			info += ";" + it->second;
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

ostream & operator << (ostream &o, VCF vcf) {

	o << vcf.Combine();
    return o;
}




