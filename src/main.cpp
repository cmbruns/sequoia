#include <iostream>
#include "AlignmentMethod.h"
#include "BioSequence.h"

using namespace std;

int main (int argc, char * const argv[]) {

	BioSequence s1 = "A";
	BioSequence s2 = "MMA";
	
	SequenceAlignment a1 = s1;
	SequenceAlignment a2 = s2;
	
	cout << align(a1, a2);
	
    return 0;
}
