#include <iostream>
using namespace std;

int main () {

	string seq1 = "aaaa";
	string seq2 = "aana";
	string seq3 = "aaNa";
	string seq4 = "aNna";
	string seq5 = "aanN";

	size_t f1 = seq1.find("nN");
	size_t f2 = seq2.find("nN");
	size_t f3 = seq3.find("nN");
	size_t f4 = seq4.find("nN");
	size_t f5 = seq5.find("nN");

	cout << seq1 << "\t" << f1 << "\n";
	cout << seq2 << "\t" << f2 << "\n";
	cout << seq3 << "\t" << f3 << "\n";
	cout << seq4 << "\t" << f4 << "\n";
	cout << seq5 << "\t" << f5 << "\n";

}
