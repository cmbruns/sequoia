// AlignmentMethod.cpp
//  routines for recursion in dynamic programming sequence alignment
//  includes class methods for classes ConservidueAlignment and AlignmentStep
//
//  $Id$
//
//  $Header$
//
//  $Log$
//  Revision 1.2  2004/05/14 15:06:47  cmbruns
//  Changed constructors
//  Added subroutines under recurrence
//  Added print statements for debugging
//  Alignment sort of works now
//
//  Revision 1.1.1.1  2004/05/11 20:26:12  cmbruns
//  Initial Mac repository for latest sequoia
//
//

#include <vector>
#include "AlignmentMethod.h"
#include "Score.h"
#include "Exceptions.h"

double BAD_SCORE = -1e50;

// ////////////////////////////// 1) Classes //////////////////////////////////////

// AlignmentStep is contained in the ConservidueAlignment class
// Each AlignmentStep represents one step in a traceback path
// of the dynamic programming table
// There is a separate AlignmentStep for each potential match, deletion, insertion,
// (and others?) state in the alignment path.
class AlignmentStep {
	friend class ConservidueAlignment;
protected:
public:
	const AlignmentStep * traceback_pointer;
	double path_score;
	const ConservidueAlignment * parent_cell;

	AlignmentStep() : 
			traceback_pointer(NULL),
			path_score(-1000),
			parent_cell(NULL)
		{}
	AlignmentStep(const ConservidueAlignment * parent) : 
		traceback_pointer(NULL),
		path_score(-1000),
		parent_cell(parent)
	{}
	
	ostream & print(ostream & os, unsigned int indent_size) const {
		string indent = "";
		for (unsigned int i=0;i<indent_size;i++) {indent += " ";}
		os << indent << "step address = " << this << endl;
		os << indent << "traceback pointer = " << traceback_pointer << endl;
		os << indent << "path_score = " << path_score << endl;
		os << indent << "parent_cell = " << parent_cell << endl;
		return os;
	}
};

// ConservidueAlignment is a single cell in a dynamic programming
// sequence alignment table.  It relates one Residue/profile in the
// first sequence, to one in the second.  For example, one ConservidueAlignment
// might represent the prospect of alignment residue VAL 123 of sperm whale
// myoglobin to residue LYS 17 if human hemoglobin.
class ConservidueAlignment {
protected:
public:
	// default constructor
	ConservidueAlignment() :
	deletion(NULL),
	insertion(NULL),
	match(NULL),
	nonmatch(NULL),
	best(NULL),
	conservidue1(NULL),
	conservidue2(NULL),
	conservidue_score(0)
	{
		// Correct self referential pointers
		deletion.parent_cell = this;
		insertion.parent_cell = this;
		match.parent_cell = this;
		nonmatch.parent_cell = this;
		best = & match;
	}

	// copy constructor
	ConservidueAlignment(const ConservidueAlignment & c2) {
		*this = c2;
	}
	
	// assignment operator
	ConservidueAlignment & operator=(const ConservidueAlignment & c2) {
		if (&c2 == this) return *this;
		
		deletion = c2.deletion;
		insertion = c2.insertion;
		match = c2.match;
		nonmatch = c2.nonmatch;
		best = c2.best;

		conservidue1 = c2.conservidue1;
		conservidue2 = c2.conservidue2;
		conservidue_score = c2.conservidue_score;
		
		// Correct self referential pointers
		deletion.parent_cell = this;
		insertion.parent_cell = this;
		match.parent_cell = this;
		nonmatch.parent_cell = this;

		if (c2.best == &c2.insertion) best = &insertion;
		if (c2.best == &c2.deletion) best = &deletion;
		if (c2.best == &c2.match) best = &match;
		if (c2.best == &c2.nonmatch) best = &nonmatch;
		
		return *this;
	}
	
	AlignmentStep deletion; // delete is a C++ keyword..., gap in sequence 1, loop in 2, "E" by Gusfield
	AlignmentStep insertion; // insertion path step(s), gap in sequence 2, loop in 1, "F" by Gusfield
	AlignmentStep match; // no gap, alignment of sequence1 to sequence 2
	AlignmentStep nonmatch; // consume both Conservidues, but do not align them (for POA)
	AlignmentStep * best; // highest scoring path step

	const Conservidue * conservidue1; // Conservidue/profile from first sequence
	const Conservidue * conservidue2; // Conservidue/profile from second sequence
	double conservidue_score; // Log odds score for aligning these two Conservidues

	void alignment_recurrence(); // Dynamic programming recurrence

	void initialize_upper_left();
	void initialize_left_column();
	void initialize_top_row();
	void compute_pair_score();
	void assign_best_match_state();
	void assign_best_deletion_state();
	void assign_best_insertion_state();
	void assign_best_best_state();
	
	const ConservidueAlignment * get_cell(const Conservidue * c1, const Conservidue * c2) const;

	ostream & print(ostream & os, unsigned int indent_size) const {
		string indent = "";
		for (unsigned int i=0;i<indent_size;i++) {indent += " ";}

		os << indent << "cell address = " << this << endl;

		os << indent << "conservidue1 = " << conservidue1 << endl;
		os << indent << "conservidue2 = " << conservidue2 << endl;
		os << indent << "conservidue_score = " << conservidue_score << endl;

		os << indent << "match step = " << endl;
		match.print(os, indent_size + 2);

		os << indent << "insertion step = " << endl;
		insertion.print(os, indent_size + 2);
		
		os << indent << "deletion step = " << endl;
		deletion.print(os, indent_size + 2);
		
		os << indent << "best step = " << best << endl;
		
		return os;
	}
};


// ////////////////////////////// 2) Functions //////////////////////////////////////

typedef vector<ConservidueAlignment> DPRow;
typedef vector<DPRow> DPMatrix;
DPMatrix * dp_matrix_pointer;
const ConservidueAlignment * 
		ConservidueAlignment::get_cell(const Conservidue * c1, 
									   const Conservidue * c2) const {
	int i = c1->array_sequence_index;
	int j = c2->array_sequence_index;
	return & (*dp_matrix_pointer)[i][j];
}

ostream & print_dp_table(const DPMatrix & mat, ostream & os, unsigned int indent_size) {
	string indent = "";
	for (unsigned int i=0;i<indent_size;i++) {indent += " ";}

	for (unsigned int i = 0; i < mat.size(); i++) {
		for (unsigned int j = 0; j < mat[i].size(); j++) {
			os << indent << "Cell(" << i << ", " << j << "):" << endl;
			mat[i][j].print(os, indent_size + 2);
			os << endl;
		}
	}

	return os;
}

// Basic sequence alignment rountine, for a single step of progressive alignment
// Should this be a member function of some class?
SequenceAlignment align(const SequenceAlignment & seq1, const SequenceAlignment & seq2) {
	SequenceAlignment answer;
	
	// Set up dynamic programming table
	int m = seq1.length();
	int n = seq2.length();

	// Initialize Dynamic programming table
	ConservidueAlignment cell;
	DPRow dp_row(n, cell);
	DPMatrix dp_table(m, dp_row);
	dp_matrix_pointer = & dp_table; // static global variable for other routines :(

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
			
			// Let the cell's methods handle the rest of the recurrence
			cell.alignment_recurrence();
		}
	}

	// For debugging:
	// print_dp_table(dp_table, cout, 2);

	// traceback
	const AlignmentStep * path_step = dp_table[m-1][n-1].best;
	cout << "Path score = " << path_step->path_score << endl;
	while (path_step != NULL) {
		const ConservidueAlignment * cell = path_step->parent_cell;
		const Conservidue * conservidue1 = cell->conservidue1;
		const Conservidue * conservidue2 = cell->conservidue2;
		int index1 = conservidue1->array_sequence_index;
		int index2 = conservidue2->array_sequence_index;
		cout << index1 << ", " << index2 << endl;
		
		path_step = path_step->traceback_pointer;
	}

	throw CBRUNS_TODO_CODE_EXCEPTION;
	
	return answer;
}

// One ConservidueAlignment is on cell of the dynamic programming table
// alignment recurrence is a member function that computes this cells contribution to the
// dynamic programming algorithm
void ConservidueAlignment::alignment_recurrence() { // Dynamic programming recurrence
	// 1) Initialization
	// Initialize begin states of each sequence (left end gap situations)

	// conservidue1->print(cout, 2);
	// cout << endl;
	// conservidue2->print(cout, 2);
	// cout << endl; // works here 
	
	// 1A Upper left corner of dp matrix, should be first cell computed
	if (conservidue1->is_initial && conservidue2->is_initial) {
		initialize_upper_left();
	}

	// 1B Possible left end gap in sequence 1
	// Top row of dp matrix
	// i == 0
	else if (conservidue1->is_initial) { // left end gap in sequence 1, is_initial is not really a Conservidue
		initialize_top_row();
	}
		
	// 1C Possible left end gap in sequence 2
	// Left column of dp matrix
	// j == 0
	else if (conservidue2->is_initial) { // left end gap in sequence 2, is_initial is not really a Conservidue
		initialize_left_column();
	}
	
	// 2) No initialization, main recursion in inner table
	// if ((!conservidue1.is_initial) && (!conservidue2.is_initial)) {
	else {
		// compute pairwise alignment score S(i,j);
		compute_pair_score();
				
		// align = V(i-1, j-1) + S(i,j) ;
		assign_best_match_state();
		
		// deletion = max (E(i, j-1), V(i, j-1) - Wg) - Ws;
		assign_best_deletion_state();
		
		// insertion = max (F(i-1, j), V(i-1, j) - Wg) - Ws;
		assign_best_insertion_state();
		
		// best = max (
		// 	& dp_cell.align
		// 	& dp_cell.deletion
		// 	& dp_cell.insertion
		// );
		assign_best_best_state();

		if (best->traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
		if (insertion.traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
		if (deletion.traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
		if (match.traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
	}				
}

// 1A Upper left corner of dp matrix, should be first cell computed
void ConservidueAlignment::initialize_upper_left() {
	SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);

	deletion.path_score = 0;
	insertion.path_score = 0;
	match.path_score = 0; // This one initialization cell does not require gaps
	
	deletion.traceback_pointer = NULL;
	insertion.traceback_pointer = NULL;
	match.traceback_pointer = NULL;
	
	if (sequence1.left_gap_factor != 0) {
		deletion.path_score += 
		sequence1.left_gap_factor * conservidue1->gap_deletion_penalty;
		deletion.path_score += 
			sequence1.left_gap_factor * conservidue2->gap_opening_penalty;
	}
	if (sequence2.left_gap_factor != 0) {
		insertion.path_score += 
		sequence2.left_gap_factor * conservidue2->gap_deletion_penalty;
		insertion.path_score += 
			sequence2.left_gap_factor * conservidue1->gap_opening_penalty;
	}
}

// 1B Possible left end gap in sequence 1
// Top row of dp matrix
// i == 0
void ConservidueAlignment::initialize_top_row() {
	// only need to initialize insertion path step
	// deletion score is implicitly handled in this initialization
	// there is no match to the start Conservidue
	// only insertion and match from inner matrix will query this cell
				
	SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	
	deletion.path_score = BAD_SCORE;
	insertion.path_score = 0;
	match.path_score = BAD_SCORE;
	
	deletion.traceback_pointer = NULL;
	match.traceback_pointer = NULL;
	
	// V(0,j) = F(0,j) = -Wg1(0) - Wo2(0) - SUM[x = 1..j](We2(x))
	
	// find best predecessor of conservidue2 (only needed for non-branched alignments)
	const ConservidueAlignment * best_previous_cell = NULL;
	double best_score = BAD_SCORE;
	vector<ConserviduePredecessor>::const_iterator prev_res;
	for (prev_res = conservidue2->predecessors.begin();
		 prev_res != conservidue2->predecessors.end();
		 prev_res ++) {
		Conservidue * previous_residue = prev_res->predecessor_conservidue;
		const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_residue);
		
		double test_score = 0;
		test_score += prev_res->transition_score; // zero for non-branched alignments
		test_score += previous_cell->insertion.path_score;
		
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			best_previous_cell = previous_cell;
		}
	}
	
	insertion.traceback_pointer = & (best_previous_cell->insertion);
	insertion.path_score = best_score;
	
	if (sequence1.left_gap_extension_factor != 0) {
		insertion.path_score += 
		sequence1.left_gap_extension_factor * conservidue2->gap_extension_penalty;
	}
}

// 1C Possible left end gap in sequence 2
// Left column of dp matrix
// j == 0
void ConservidueAlignment::initialize_left_column() {
	// only need to initialize deletion path step
	// insertion score is implicitly handled in this initialization
	// there is no match to the start Conservidue
	// only deletion and match from inner matrix will query this cell
	
	SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);

	insertion.path_score = BAD_SCORE;
	deletion.path_score = 0;
	match.path_score = BAD_SCORE;
	
	insertion.traceback_pointer = NULL;
	match.traceback_pointer = NULL;
	
	// V(i,0) = E(i,0) = -Wg2(0) - Wo1(0) - SUM[x = 1..i](We1(x))
	// TODO 
	
	// find best predecessor of conservidue1 (only needed for non-branched alignments)
	const ConservidueAlignment * best_previous_cell = NULL;
	double best_score = BAD_SCORE;
	vector<ConserviduePredecessor>::const_iterator prev_res;
	for (prev_res = conservidue1->predecessors.begin();
		 prev_res != conservidue1->predecessors.end();
		 prev_res ++) {
		Conservidue * previous_residue = prev_res->predecessor_conservidue;
		const ConservidueAlignment * previous_cell = get_cell(previous_residue, conservidue2);
		
		double test_score = 0;
		test_score += prev_res->transition_score; // zero for non-branched alignments
		test_score += previous_cell->deletion.path_score;
		
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			best_previous_cell = previous_cell;
		}
	}
	
	deletion.traceback_pointer = & (best_previous_cell->deletion);
	deletion.path_score = best_score;
	
	if (sequence2.left_gap_extension_factor != 0) {
		deletion.path_score += 
		sequence2.left_gap_extension_factor * conservidue1->gap_extension_penalty;
	}
}

// compute pairwise alignment score S(i,j);
void ConservidueAlignment::compute_pair_score() {
	conservidue_score = conservidue_pair_score(*conservidue1, *conservidue2);
}

// align = V(i-1, j-1) + S(i,j) ;
// G(i,j) = max{
// 	E(i-1, j-1) - Wc2(j-1), 
// 	F(i-1, j-1) - Wc1(i-1), 
// 	G(i-1, j-1)
// } + S(i,j) 
void ConservidueAlignment::assign_best_match_state() {
	double best_score = BAD_SCORE;
	AlignmentStep * subject_step = & match;

	// examine all predecessor pairs
	vector<ConserviduePredecessor>::const_iterator prev_res1;
	vector<ConserviduePredecessor>::const_iterator prev_res2;
	for (prev_res1 = conservidue1->predecessors.begin();
		 prev_res1 != conservidue1->predecessors.end();
		 prev_res1 ++) {
		
		const ConserviduePredecessor & cp1 = * prev_res1;
		const Conservidue * previous_conservidue1 = cp1.predecessor_conservidue;
		float transition_score1 = cp1.transition_score;
		
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {

			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = cp2.predecessor_conservidue;
			float transition_score2 = cp2.transition_score;
			
			const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, previous_conservidue2);
			double test_score;
			const AlignmentStep * test_step;

			// Beware of special case of closing an initial end gap
			float end_gap_factor1 = 1.0;
			float end_gap_factor2 = 1.0;
			if (previous_conservidue2->is_initial)
				end_gap_factor2 = conservidue2->parent_alignment->left_gap_factor;
			if (previous_conservidue1->is_initial)
				end_gap_factor1 = conservidue1->parent_alignment->left_gap_factor;
			
			// A) match to match
			// G(i-1, j-1) + S(i,j) 
			test_step = & previous_cell->match;
			test_score = test_step->path_score + 
				transition_score1 +
				transition_score2 +
				conservidue_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}

			// B) match to deletion
			// 	E(i-1, j-1) - Wc2(j-1) + S(i,j)
			test_step = & previous_cell->deletion;
			test_score = test_step->path_score + 
				transition_score1 +
				transition_score2 +
				previous_conservidue2->gap_closing_penalty * end_gap_factor2 +
				conservidue_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}

			// C) match to insertion
			// 	F(i-1, j-1) - Wc1(i-1) + S(i,j)
			test_step = & previous_cell->insertion;
			test_score = test_step->path_score + 
				transition_score1 +
				transition_score2 +
				previous_conservidue1->gap_closing_penalty * end_gap_factor1 +
				conservidue_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		}
	}
	if (subject_step->traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
}

// E(i,j) = max {
// 	E(i, j-1), 
// 	G(i, j-1) - Wg1(i) - Wo2(j-1), 
// 	F(i, j-1) - Wg1(i) - Wo2(j-1) - Wc1(i) (?)
// } - We2(j)
void ConservidueAlignment::assign_best_deletion_state() {
	double best_score = BAD_SCORE;
	AlignmentStep * subject_step = & deletion;
	
	float end_gap_factor = 1.0; // In case end gaps are free
	float end_gap_extension_factor = 1.0; // In case end gaps are free
	if (conservidue1->is_final) {
		end_gap_factor = conservidue1->parent_alignment->right_gap_factor;
		end_gap_extension_factor = conservidue1->parent_alignment->right_gap_extension_factor;
	}
	// Perhaps a final end gap closing penalty is needed
	float final_gap_close = 0.0;
	if (conservidue1->is_final && conservidue2->is_final) final_gap_close = conservidue2->gap_closing_penalty * end_gap_factor;
	
	// examine all predecessors
	vector<ConserviduePredecessor>::const_iterator prev_res2;
	for (prev_res2 = conservidue2->predecessors.begin();
		prev_res2 != conservidue2->predecessors.end();
		prev_res2 ++) {
			
		const ConserviduePredecessor & cp2 = * prev_res2;
		const Conservidue * previous_conservidue2 = cp2.predecessor_conservidue;
		float transition_score2 = cp2.transition_score;
			
		const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_conservidue2);
		double test_score;
		const AlignmentStep * test_step;
			
		// A) deletion to deletion
		// 	E(i, j-1) - We2(j)
		test_step = & previous_cell->deletion;
		test_score = test_step->path_score + 
			transition_score2 +
			final_gap_close +
			conservidue2->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
		
		// B) deletion to match
		// 	G(i, j-1) - Wg1(i) - Wo2(j-1) - We2(j)
		test_step = & previous_cell->match;
		test_score = test_step->path_score + 
			transition_score2 +
			final_gap_close +
			conservidue1->gap_deletion_penalty * end_gap_factor +
			previous_conservidue2->gap_opening_penalty * end_gap_factor +
			conservidue2->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
		
		// C) deletion to insertion
		// 	F(i, j-1) - Wg1(i) - Wo2(j-1) - Wc1(i) - We2(j)
		test_step = & previous_cell->insertion;
		test_score = test_step->path_score + 
			transition_score2 +
			final_gap_close +
			conservidue1->gap_deletion_penalty * end_gap_factor +
			previous_conservidue2->gap_opening_penalty * end_gap_factor +
			conservidue1->gap_closing_penalty + // Not part of any end gap!
			conservidue2->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
	}		
	if (subject_step->traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
}

// F(i,j) = max {
// 	F(i-1, j), 
// 	G(i-1, j) - Wg2(j) -Wo1(i-1), 
// 	E(i-1, j) - Wg2(j) -Wo1(i-1) - Wc2(j) (?)
// } - We1(i)
void ConservidueAlignment::assign_best_insertion_state() {
	double best_score = BAD_SCORE;
	AlignmentStep * subject_step = & insertion;

	float end_gap_factor = 1.0; // In case end gaps are free
	float end_gap_extension_factor = 1.0; // In case end gaps are free
	if (conservidue2->is_final) {
		end_gap_factor = conservidue2->parent_alignment->right_gap_factor;
		end_gap_extension_factor = conservidue2->parent_alignment->right_gap_extension_factor;
	}
	// Perhaps a final end gap closing penalty is needed
	float final_gap_close = 0.0;
	if (conservidue1->is_final && conservidue2->is_final) final_gap_close = conservidue1->gap_closing_penalty * end_gap_factor;
	
	// examine all predecessors
	vector<ConserviduePredecessor>::const_iterator prev_res1;
	for (prev_res1 = conservidue1->predecessors.begin();
		 prev_res1 != conservidue1->predecessors.end();
		 prev_res1 ++) {
		
		const ConserviduePredecessor & cp1 = * prev_res1;
		const Conservidue * previous_conservidue1 = cp1.predecessor_conservidue;
		float transition_score1 = cp1.transition_score;
		
		const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, conservidue2);
		double test_score;
		const AlignmentStep * test_step;
		
		// A) insertion to insertion
		// 	F(i-1, j) - We1(j)
		test_step = & previous_cell->insertion;
		test_score = test_step->path_score + 
			transition_score1 +
			final_gap_close +
			conservidue1->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
		
		// B) insertion to match
		// 	G(i-1, j) - Wg2(j) -Wo1(i-1) - We1(j)
		test_step = & previous_cell->match;
		test_score = test_step->path_score + 
			transition_score1 +
			final_gap_close +
			conservidue2->gap_deletion_penalty * end_gap_factor +
			previous_conservidue1->gap_opening_penalty * end_gap_factor +
			conservidue1->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
		
		// C) insertion to deletion
		// 	E(i-1, j) - Wg2(j) -Wo1(i-1) - Wc2(j) - We1(i)
		test_step = & previous_cell->deletion;
		test_score = test_step->path_score + 
			transition_score1 +
			final_gap_close +
			conservidue2->gap_deletion_penalty * end_gap_factor +
			previous_conservidue1->gap_opening_penalty * end_gap_factor +
			conservidue2->gap_closing_penalty + // This one cannot be part of end gap!
 			conservidue1->gap_extension_penalty * end_gap_extension_factor;
		if ((test_score > best_score) || (best_score == BAD_SCORE)) {
			best_score = test_score;
			subject_step->path_score = test_score;
			subject_step->traceback_pointer = test_step;
		}
	}		
	if (subject_step->traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
}

// best = max (
// 	& dp_cell.align
// 	& dp_cell.deletion
// 	& dp_cell.insertion
// );
void ConservidueAlignment::assign_best_best_state() {
	best = & match;
	if (insertion.path_score > best->path_score) best = & insertion;
	if (deletion.path_score > best->path_score) best = & deletion;
	if (best->traceback_pointer == NULL) throw CBRUNS_DEBUGGING_EXCEPTION;
}


