#include <iostream>
#include "AlignmentMethod.h"
#include "BioSequence.h"

using namespace std;

int main (int argc, char * const argv[]) {

	BioSequence s1 = "ADGH";
	BioSequence s2 = "GHKLADGHL";
	cout << s1 << endl;
	
	SequenceAlignment a1 = s1;
	SequenceAlignment a2 = s2;
	
	cout << "Test of alignment structures" << endl;
	cout << a1 << endl;
	cout << a2 << endl;
	cout << align(a1, a2); // TODO
	
    return 0;
}
