#include <iostream>
#include <fstream>
#include "AlignmentMethod.h"
#include "BioSequence.h"
#include "Exceptions.h"

using namespace std;

int main (int argc, char * const argv[]) {

	// Read first sequence file
	const char * infile1_name = "seq1.fasta";
	ifstream infile1(infile1_name);
	if (infile1 == 0) {
		cerr << "*** ERROR *** : Unable to open file " << infile1_name << endl;
		return NO_SUCH_FILE_EXCEPTION;
	}
	SequenceAlignment a1;
	infile1 >> a1;
	
	cout << a1;
	
	// Read second sequence file
	const char * infile2_name = "seq2.fasta";
	ifstream infile2(infile2_name);
	if (infile2 == 0) {
		cerr << "*** ERROR *** : Unable to open file " << infile2_name << endl;
		return NO_SUCH_FILE_EXCEPTION;
	}
	SequenceAlignment a2;
	infile2 >> a2;
	
	cout << a2;

	// Write aligned sequence file
	const char * outfile_name = "align.fasta";
	ofstream outfile (outfile_name);
	if (outfile == 0) {
		cerr << "*** ERROR *** : Unable to open file " << outfile_name << endl;
		return NO_SUCH_FILE_EXCEPTION;
	}
	SequenceAlignment a3 = a2.align(a1);
	outfile << a3; // The heavy lifting is in this line

	cout << a3;
	
    return 0;
}
