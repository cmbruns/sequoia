#include <iostream>
#include <fstream>
#include "AlignmentMethod.h"
#include "BioSequence.h"
#include "Exceptions.h"

using namespace std;

int main (int argc, char * const argv[]) {

	BioSequence s1 = "AAWWAACCAA";
	BioSequence s2 = "WWCC";
	BioSequence s3 = "HHHMITNLERVGHPKERCFTF";
	SequenceAlignment align1 = s1;	
	SequenceAlignment align2 = s2;
	SequenceAlignment align3 = s3;
	align1.set_weight(0.20);
	align2.set_weight(1.65);
	align3.set_weight(0.96);
	
	SequenceAlignment true_alignment;
	true_alignment.load_fasta_file("fnr_pdr_true.fasta");
	// cout << true_alignment << endl;

	// Read first sequence file
	SequenceAlignment a1;
	a1.load_fasta_file("fnr.fasta");
	a1.set_weight(1.0);

	// Read second sequence file
	SequenceAlignment a2;
	a2.load_fasta_file("pdr.fasta");
	a2.set_weight(1.0);
	
	SequenceAlignment a3;
	a3.load_fasta_file("b5.fasta");
	a3.set_weight(1.0);
	
	AlignmentScore cumulative_alignment_score;
	
	SequenceAlignment alignment;
	alignment = a1.align(a2);
	AlignmentScore SOP = alignment.sum_of_pairs_score();
	// align1.print_debug(cout);
	// align2.print_debug(cout);
	// alignment.print_debug(cout);

	cout << alignment << endl;
	cout << "Sum of pairs score = " << SOP << endl;
	SOP.print_details(cout);
	cumulative_alignment_score += alignment.alignment_score();
	cout << "Cumulative alignment score = " << cumulative_alignment_score << endl;
	cumulative_alignment_score.print_details(cout);
	cout << "Discrepancy = " << SOP - cumulative_alignment_score << endl;
	cout << endl << endl;
	
	// alignment.print_debug(cout);
	// align3.print_debug(cout);
	alignment = a3.align(alignment); // Add third sequence to alignment
	SOP = alignment.sum_of_pairs_score();
	// alignment.print_debug(cout);

	cout << alignment << endl;
	cout << "Sum of pairs score = " << SOP << endl;
	SOP.print_details(cout);
	cumulative_alignment_score += alignment.alignment_score();
	cout << endl << endl;
	cout << "Cumulative alignment score = " << cumulative_alignment_score << endl;
	cout << "Discrepancy = " << SOP - cumulative_alignment_score << endl;
	cumulative_alignment_score.print_details(cout);
	cout << endl << endl;
	
	// cout << "Sensitivity = " << alignment.report_accuracy(true_alignment) << endl;

    return 0;
}
