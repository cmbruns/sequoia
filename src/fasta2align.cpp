/* Copyright (c) 2005 Christopher M. Bruns
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// $Id$
// $Header$
// $Log$
// Revision 1.6  2004/06/14 16:59:18  cmbruns
// Reworked yet again for different testing of sum of pairs score
//
// Revision 1.5  2004/06/04 19:32:24  cmbruns
// Updated GPL header
//
// Changed for initial testing of structure accessible surface area computation.  This has been extended in residue_area.cpp
//


#include <iostream>
#include <fstream>
#include "SequenceAlignment.h"
#include "BioSequence.h"
#include "Exceptions.h"

using namespace std;

int main (int argc, char * const argv[]) {

	SequenceAlignment b1,b2,b3;
	
	
	// SequenceAlignment true_alignment;
	// true_alignment.load_fasta_file("fnr_pdr_true.fasta");
	// cout << true_alignment << endl;
		
	b1.load_pdb_file("1FND.pdb");
	// b1.set_weight(1.0);
	
	b2.load_fasta_file("pdr.fasta");
	// b2.set_weight(0.7);

	b3.load_fasta_file("b5.fasta");
	// b3.set_weight(0.7);

	// b1 = "AAACCCWWW";
	// b3 = "CCCWWW";
	// b2 = "CCCWWW";
	
	SequenceAlignment test_alignment = b1.align(b2);
	cout << test_alignment;
	
	// cout << "Accuracy = " << test_alignment.report_accuracy(true_alignment) << endl;

	AlignmentScore alignment_score = test_alignment.alignment_score();
	AlignmentScore sum_of_pairs_score = test_alignment.delta_sum_of_pairs_score(b1.sequence().size());
	AlignmentScore discrepancy_score = alignment_score - sum_of_pairs_score;
	
	cout << endl << endl;
	cout << "Alignment score:" << endl;
	alignment_score.print_details(cout);
	cout << endl << endl;
	cout << "Sum of pairs score:" << endl;
	sum_of_pairs_score.print_details(cout);
	cout << endl << endl;
	cout << "Discrepancy:" << endl;
	discrepancy_score.print_details(cout);
	cout << endl << endl;
	

	
	SequenceAlignment test_alignment2 = b3.align(test_alignment);
	cout << test_alignment2;
	
	// cout << "Accuracy = " << test_alignment2.report_accuracy(true_alignment) << endl;
	
	alignment_score = test_alignment2.alignment_score();
	sum_of_pairs_score = test_alignment2.delta_sum_of_pairs_score(b3.sequence().size());
	discrepancy_score = alignment_score - sum_of_pairs_score;
	
	cout << endl << endl;
	cout << "Alignment score:" << endl;
	alignment_score.print_details(cout);
	cout << endl << endl;
	cout << "Sum of pairs score:" << endl;
	sum_of_pairs_score.print_details(cout);
	cout << endl << endl;
	cout << "Discrepancy:" << endl;
	discrepancy_score.print_details(cout);

	// cout << endl << endl << "first alignment" << endl;
	// test_alignment.print_debug();

	// cout << endl << endl << "final alignment" << endl;
	// test_alignment2.print_debug();
	
	return 0;
}
