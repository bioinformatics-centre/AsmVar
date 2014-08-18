#include <iostream>
#include <sstream>

using namespace std;

void test(int n, int m) {
	
	cout << "Number : " << n << "\tRaw: " << m << endl;
}

int main() {

	double dbl = 3.1415926;

	ostringstream strs;
	strs << dbl;
	string str = strs.str();
	cout << str << "\n";
	for (int a = 0; a < 20; ++a ) test( a > 10 ? a : 8, a );
}
