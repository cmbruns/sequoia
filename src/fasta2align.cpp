// This file is part of the Sequoia package for macromolecular 
//  sequence/structure analysis
// Copyright (C) 2004  Christopher M. Bruns, Ph.D.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// 
// See the accompanying file 'LICENSE' for details
// 
// To contact the author, write to cmbruns@comcast.net or bruns@scripps.edu
// In publications please cite: Bruns et al (1999), J.Mol.Biol. 288:427-439
// Please submit bug reports at http://bruns.homeip.net/bugzilla/
// 
// To obtain a non-GPL version of this program, see http://bruns.homeip.net/sequoia.html
//

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
