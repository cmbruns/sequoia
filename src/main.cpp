#include <iostream>
#include "AlignmentMethod.h"
#include "BioSequence.h"

using namespace std;

int main (int argc, char * const argv[]) {

	BioSequence s1 = "A";
	BioSequence s2 = "A";
	
	cout << s1 << endl;
	cout << s2 << endl;
	
	SequenceAlignment a1 = s1;
	SequenceAlignment a2 = s2;
	
	cout << a1 << endl;
	cout << a2 << endl;
	
	cout << a1.align(a2);
	
    return 0;
}
