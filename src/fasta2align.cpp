#include <iostream>
#include <fstream>
#include "AlignmentMethod.h"
#include "BioSequence.h"

using namespace std;

int main (int argc, char * const argv[]) {
	
	SequenceAlignment a1;
	SequenceAlignment a2;
	
	ifstream file1("seq1.fasta");
	ifstream file2("seq2.fasta");
	
	file1 >> a1;
	file2 >> a2;
	
	ofstream outfile ("align.fasta");
	
	outfile << align(a1, a2);
	
    return 0;
}
