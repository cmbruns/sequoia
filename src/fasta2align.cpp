#include <iostream>
#include <fstream>
#include "AlignmentMethod.h"
#include "BioSequence.h"
#include "Exceptions.h"

using namespace std;

int main (int argc, char * const argv[]) {

	SequenceAlignment true_alignment;
	ifstream true_infile("fnr_pdr_true.fasta");
	true_infile >> true_alignment;
	true_infile.close();
	cout << true_alignment << endl;

	// Read first sequence file
	const char * infile1_name = "seq1.fasta";
	ifstream infile1(infile1_name);
	if (infile1 == 0) {
		cerr << "*** ERROR *** : Unable to open file " << infile1_name << endl;
		return NO_SUCH_FILE_EXCEPTION;
	}
	SequenceAlignment a1;
	infile1 >> a1;

	// Read second sequence file
	const char * infile2_name = "seq2.fasta";
	ifstream infile2(infile2_name);
	if (infile2 == 0) {
		cerr << "*** ERROR *** : Unable to open file " << infile2_name << endl;
		return NO_SUCH_FILE_EXCEPTION;
	}
	SequenceAlignment a2;
	infile2 >> a2;
	infile2.close();
	
	SequenceAlignment a3;
	a3 = a1.align(a2, ALIGN_LOCAL);
	cout << a3 << endl;
	cout << "Sensitivity = " << a3.report_accuracy(true_alignment) << endl;

    return 0;
}
