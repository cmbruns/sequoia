// AlignmentMethod.cpp
//  routines for recursion in dynamic programming sequence alignment
//  includes class methods for classes ConservidueAlignment and AlignmentStep
//
//  $Id$
//  $Header$
//  $Log$
//  Revision 1.1  2004/05/11 20:26:12  cmbruns
//  Initial revision
//
//

#include "AlignmentMethod.h"
#include "Score.h"

// ////////////////////////////// 1) Classes //////////////////////////////////////

// AlignmentStep is contained in the ConservidueAlignment class
// Each AlignmentStep represents one step in a traceback path
// of the dynamic programming table
// There is a separate AlignmentStep for each potential match, insert, delete,
// (and others?) state in the alignment path.
class AlignmentStep {
	AlignmentStep * traceback_pointer;
	double path_score;
	bool is_aligned;  // What is this for?
					  // ConservidueAlignment * parent; // Is this needed?
};

// ConservidueAlignment is a single cell in a dynamic programming
// sequence alignment table.  It relates one Residue/profile in the
// first sequence, to one in the second.  For example, one ConservidueAlignment
// might represent the prospect of alignment residue VAL 123 of sperm whale
// myoglobin to residue LYS 17 if human hemoglobin.
class ConservidueAlignment {
public:
	AlignmentStep insertion; // insert path step(s), gap in sequence 1
	AlignmentStep deletion; // delete is a C++ keyword..., gap in sequence 2
	AlignmentStep match; // no gap, alignment of sequence1 to sequence 2
	AlignmentStep nonalign; // consume both Conservidues, but do not align them (for POA)
	AlignmentStep * best; // highest scoring path step
	const Conservidue * conservidue1; // Conservidue/profile from first sequence
	const Conservidue * conservidue2; // Conservidue/profile from second sequence
	double conservidue_score; // Log odds score for aligning these two Conservidues
	void alignment_recurrence(); // Dynamic programming recurrence
};


// ////////////////////////////// 2) Functions //////////////////////////////////////

// Basic sequence alignment rountine, for a single step of progressive alignment
// Should this be a member function of some class?
SequenceAlignment align(const SequenceAlignment & seq1, const SequenceAlignment & seq2) {
	SequenceAlignment answer;
	
	// Set up dynamic programming table
	int m = seq1.length();
	int n = seq2.length();
	ConservidueAlignment dp_table[m][n];

	// Visit each cell of the table and run alignment recurrence
	int i, j;
	const Conservidue * residue1;
	const Conservidue * residue2;
	for (i = 0; i < m; ++i) {
		residue1 = & seq1[i];
		for (j = 0; j < n; ++j) {
			residue2 = & seq2[j];

			// Initialize one cell of the table
			ConservidueAlignment & cell = dp_table[i][j];
			cell.conservidue1 = residue1;
			cell.conservidue2 = residue2;
			cell.conservidue_score = conservidue_score(*residue1, *residue2);
			
			// Let the cell's methods handle the rest of the recurrence
			cell.alignment_recurrence();
		}
	}

	// TODO traceback
	throw 100;
	
	return answer;
}

void alignment_recurrence() { // Dynamic programming recurrence
	// Initialize
	if (conservidue2.is_begin) // left end gap in sequence 2, is_begin is not really a Conservidue
// only need to initialize insert path step
// delete score is implicitly handled in this initialization
// there is no match to the start Conservidue
// only insert and match from inner matrix will query this cell
		
		double & score = dp_cell.insert.path_score;
	if (left_extension_penalty_factor_2 == 0)
		score = 0;
	else 
		score = left_extension_penalty_factor_2 * 
			extension_penalty(conservidue1.begin_distance)
			score += conservidue2.gap_close_penalty // even with end gap free.  If this should be zero, set it to zero at the sequence initialization stage
			score += conservidue1.gap_close_penalty // even with end gaps (otherwise) free
			
			dp_cell.best = & dp_cell.insert // the only valid state
			
			if (conservidue1.is_begin) // left end gap in sequence 1, is_begin is not really a Conservidue
// only need to initialize delete path step
// insert score is implicitly handled in this initialization
// there is no match to the start Conservidue
// only delete and match from inner matrix will query this cell
				
				double & score = dp_cell.delete.path_score;
	if (left_extension_penalty_factor_1 == 0)
		score = 0;
	else 
		score = left_extension_penalty_factor_1 * 
			extension_penalty(conservidue2.begin_distance)
			score += conservidue1.gap_close_penalty // even with end gap free.  If this should be zero, set it to zero at the sequence initialization stage
			score += conservidue2.gap_close_penalty // even with end gaps (otherwise) free
			
			dp_cell.best = & dp_cell.delete // the only valid state
			
// No initialization, main recursion in inner table
			if ((!conservidue1.is_begin) && (!conservidue2.is_begin)) {
				compute pairwise alignment score S(i,j)
				
				dp_cell.align = V(i-1, j-1) + S(i,j) 
				dp_cell.insert = max {E(i, j-1), V(i, j-1) - Wg} - Ws
				dp_cell.delete = max {F(i-1, j), V(i-1, j) - Wg} - Ws
				dp_cell.best = max {
					& dp_cell.align
					& dp_cell.insert
					& dp_cell.delete
				}
				
	throw 100; // TODO
}

