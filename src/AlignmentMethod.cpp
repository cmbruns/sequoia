// AlignmentMethod.cpp
//  routines for recursion in dynamic programming sequence alignment
//  includes class methods for classes ConservidueAlignment and AlignmentStep
//
//  $Id$
//
//  $Header$
//
//  $Log$
//  Revision 1.3  2004/05/24 14:52:06  cmbruns
//  For all insertion, deletion, and nonmatch states, iterate over each segment of piecewise linear gap model
//  ALERT: Default gap segment count is a kludge and needs to be revisited TODO
//  Added is_match member to ConservidueAlignment
//  Moved ConservidueAlignment constructors to end of block
//  Added parent_table member to ConservidueAlignment: no more global variable for finding dynamic programming table
//  Implemented local alignment option
//    add optional granularity argument to most alignment subroutines
//  Completed traceback portion of align subroutine
//
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

#define DEFAULT_GAP_SEGMENT_COUNT 5
double BAD_SCORE = -1e50;

// ////////////////////////////// 1) Classes //////////////////////////////////////

typedef vector<ConservidueAlignment> DPRow;
typedef vector<DPRow> DPMatrix;

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
	bool is_match; // true if this step actually aligns two residues, i.e. match state

	AlignmentStep() : 
			traceback_pointer(NULL),
			path_score(-1000),
			parent_cell(NULL),
			is_match(false)
		{}
	AlignmentStep(const ConservidueAlignment * parent) : 
		traceback_pointer(NULL),
		path_score(-1000),
		parent_cell(parent),
		is_match(false)
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
	vector<AlignmentStep> deletion; // delete is a C++ keyword..., gap in sequence 1, loop in 2, "E" by Gusfield
	vector<AlignmentStep> insertion; // insertion path step(s), gap in sequence 2, loop in 1, "F" by Gusfield
	vector<AlignmentStep> nonmatch; // consume both Conservidues, but do not align them (for POA)
	AlignmentStep match; // no gap, alignment of sequence1 to sequence 2
	AlignmentStep * best; // highest scoring path step
	
	DPMatrix * parent_table;
	const Conservidue * conservidue1; // Conservidue/profile from first sequence
	const Conservidue * conservidue2; // Conservidue/profile from second sequence
	double conservidue_score; // Log odds score for aligning these two Conservidues
	
	void alignment_recurrence(AlignmentGranularity granularity = ALIGN_GLOBAL); // Dynamic programming recurrence
	
	void initialize_upper_left(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void initialize_left_column(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void initialize_top_row(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void compute_pair_score();
	void assign_best_match_state(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void assign_best_deletion_state(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void assign_best_insertion_state(AlignmentGranularity granularity = ALIGN_GLOBAL);
	void assign_best_best_state(AlignmentGranularity granularity = ALIGN_GLOBAL);
	
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
		
		for (unsigned int i = 0; i < insertion.size(); i++) {
			os << indent << "insertion step " << i << " = " << endl;
			insertion[i].print(os, indent_size + 2);
		}

		for (unsigned int i = 0; i < deletion.size(); i++) {
			os << indent << "deletion step " << i << " = " << endl;
			deletion[i].print(os, indent_size + 2);
		}
		
		os << indent << "best step = " << best << endl;
		
		return os;
	}

	// default constructor
	ConservidueAlignment() :
	deletion(DEFAULT_GAP_SEGMENT_COUNT),
	insertion(DEFAULT_GAP_SEGMENT_COUNT),
	nonmatch(DEFAULT_GAP_SEGMENT_COUNT),
	match(NULL),
	best(NULL),
	conservidue1(NULL),
	conservidue2(NULL),
	conservidue_score(0)
	{
		// Correct self referential pointers
		for (unsigned int i = 0; i < deletion.size(); i ++) { 
			deletion[i].parent_cell = this;
			insertion[i].parent_cell = this;
			nonmatch[i].parent_cell = this;

			deletion[i].is_match = false;
			insertion[i].is_match = false;
			nonmatch[i].is_match = false;
		}
		match.parent_cell = this;
		match.is_match = true;

		best = & match; // by default
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
		for (unsigned int i = 0; i < deletion.size(); i ++) { 
			deletion[i].parent_cell = this;
			insertion[i].parent_cell = this;
			nonmatch[i].parent_cell = this;
		}
		match.parent_cell = this;
		
		for (unsigned int i = 0; i < deletion.size(); i ++) { 
			if (c2.best == &c2.insertion[i]) best = &insertion[i];
			if (c2.best == &c2.deletion[i]) best = &deletion[i];
			if (c2.best == &c2.nonmatch[i]) best = &nonmatch[i];
		}
		if (c2.best == &c2.match) best = &match;
		
		return *this;
	}

};


// ////////////////////////////// 2) Functions //////////////////////////////////////

// DPMatrix * dp_matrix_pointer;
const ConservidueAlignment * 
  ConservidueAlignment::get_cell(const Conservidue * c1, 
								 const Conservidue * c2) const {
 	int i = c1->array_sequence_index;
 	int j = c2->array_sequence_index;
 	return & (*parent_table)[i][j];
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
SequenceAlignment SequenceAlignment::align
		(const SequenceAlignment & seq2,
		 AlignmentGranularity granularity) const 
{
	const SequenceAlignment & seq1 = *this;

	// Set up dynamic programming table
	unsigned int m = seq1.length();
	unsigned int n = seq2.length();

	// Initialize Dynamic programming table
	ConservidueAlignment cell;
	DPRow dp_row(n, cell);
	DPMatrix dp_table(m, dp_row);
	// dp_matrix_pointer = & dp_table; // static global variable for other routines :(

	// Keep track of highest scoring spot in table, for local alignment
	AlignmentStep * maximal_segment_pair_end = NULL;
	double maximal_segment_pair_score = BAD_SCORE;

	// Visit each cell of the table and run alignment recurrence
	unsigned int i, j;
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
			cell.parent_table = & dp_table;
			
			// Let the cell's methods handle the rest of the recurrence
			cell.alignment_recurrence(granularity);

			// In support of local Smith-Waterman alignment
			if (cell.match.path_score > maximal_segment_pair_score) {
				maximal_segment_pair_score = cell.match.path_score;
				maximal_segment_pair_end = & cell.match;
			}
		}
	}

	// For debugging:
	// print_dp_table(dp_table, cout, 2);

	// traceback
	const AlignmentStep * path_step = NULL;

	if (granularity == ALIGN_LOCAL) // Local Smith-Waterman alignment
		path_step = maximal_segment_pair_end;
	else if ((m > 0) && (n > 0)) path_step = dp_table[m-1][n-1].best;
	
	double path_score = 0;
	if (path_step != NULL) path_score = path_step->path_score;
	// cout << "Alignment score = " << path_score << endl;
	
	// Traceback is in the opposite direction of convenient alignment construction,
	// so store coordinates for reverse iteration
	vector<long> alignment1_matches;
	vector<long> alignment2_matches;
	bool dump_traceback_coordinates = false; // for debugging
	while (path_step != NULL) {
		int index1 = path_step->parent_cell->conservidue1->array_sequence_index;
		int index2 = path_step->parent_cell->conservidue2->array_sequence_index;
		if (dump_traceback_coordinates) cout << index1 << ", " << index2;
		if (path_step->is_match) {
			if (dump_traceback_coordinates) cout << "*";
			alignment1_matches.push_back(index1);
			alignment2_matches.push_back(index2);
		}
		if (dump_traceback_coordinates) cout << endl;
		path_step = path_step->traceback_pointer;
	}

	// Start creating final sequence alignment structure
	SequenceAlignment answer = seq1;
	// Store the alignment score while we have it handy
	answer.pair_alignment_score = path_score;
	
	// 1) create sequences in final alignment
	answer.sequences.clear();
	for (unsigned int sequence_number = 0; 
		 sequence_number < seq1.sequences.size(); 
		 sequence_number ++) {
		answer.add_sequence(seq1.sequences[sequence_number]);
	}
	for (unsigned int seq2_sequence_number = 0; 
		 seq2_sequence_number < seq2.sequences.size(); 
		 seq2_sequence_number ++) {
		answer.add_sequence(seq2.sequences[seq2_sequence_number]);
	}
	
	// 2) TODO - update residues, if necessary
	// 3) TODO - update conservidues

	// Now step through matches in forward direction
	// TODO - update sequence_residues mapping in new conservidues
	int sequence2_index_offset = seq1.sequences.size();
	answer.conservidues.clear();
	Conservidue begin_conservidue;
	answer.add_conservidue(begin_conservidue);
	
	unsigned long alignment1_position = 0;
	unsigned long alignment2_position = 0;
	for (long match = alignment1_matches.size() - 1; match >= 0; match--) {
		unsigned long match1 = alignment1_matches[match];
		unsigned long match2 = alignment2_matches[match];
		// Insert gaps in sequence 2
		while (alignment1_position < match1) {
			if (! seq1[alignment1_position].is_initial)
				answer.add_conservidue(seq1[alignment1_position], 0);
			alignment1_position ++;
		}
		// Insert gaps in sequence 1
		while (alignment2_position < match2) {
			if (! seq2[alignment2_position].is_initial)
				answer.add_conservidue(seq2[alignment2_position], sequence2_index_offset);
			alignment2_position ++;
		}
		// insert match
		if ((! seq2[alignment2_position].is_initial) && (! seq1[alignment1_position].is_initial))
			answer.add_conservidue(seq1[alignment1_position].combine_conservidues(seq2[alignment2_position]));

		alignment1_position ++;
		alignment2_position ++;
	}
	// extend gaps to end
	// Insert gaps in sequence 2
	while (alignment1_position < m) {
		if (seq1[alignment1_position].is_initial) {
			alignment1_position ++;
			continue;
		}
		answer.add_conservidue(seq1[alignment1_position], 0);
		alignment1_position ++;
	}
	// Insert gaps in sequence 1
	while (alignment2_position < n) {
		if (seq2[alignment2_position].is_initial) {
			alignment2_position ++;
			continue;
		}
		answer.add_conservidue(seq2[alignment2_position], sequence2_index_offset);
		alignment2_position ++;
	}
	
	// throw CBRUNS_TODO_CODE_EXCEPTION;
	// TODO - reweight sequences in output?

	// answer.print_debug(); // for debugging
	
	return answer;
}

// One ConservidueAlignment is on cell of the dynamic programming table
// alignment recurrence is a member function that computes this cells contribution to the
// dynamic programming algorithm
void ConservidueAlignment::alignment_recurrence(AlignmentGranularity granularity) { // Dynamic programming recurrence
	// 1) Initialization
	// Initialize begin states of each sequence (left end gap situations)

	// conservidue1->print(cout, 2);
	// cout << endl;
	// conservidue2->print(cout, 2);
	// cout << endl; // works here 
	
	// 1A Upper left corner of dp matrix, should be first cell computed
	if (conservidue1->is_initial && conservidue2->is_initial) {
		initialize_upper_left(granularity);
	}

	// 1B Possible left end gap in sequence 1
	// Top row of dp matrix
	// i == 0
	else if (conservidue1->is_initial) { // left end gap in sequence 1, is_initial is not really a Conservidue
		initialize_top_row(granularity);
	}
		
	// 1C Possible left end gap in sequence 2
	// Left column of dp matrix
	// j == 0
	else if (conservidue2->is_initial) { // left end gap in sequence 2, is_initial is not really a Conservidue
		initialize_left_column(granularity);
	}
	
	// 2) No initialization, main recursion in inner table
	// if ((!conservidue1.is_initial) && (!conservidue2.is_initial)) {
	else {
		// compute pairwise alignment score S(i,j);
		compute_pair_score();
				
		// align = V(i-1, j-1) + S(i,j) ;
		assign_best_match_state(granularity);
		
		// deletion = max (E(i, j-1), V(i, j-1) - Wg) - Ws;
		assign_best_deletion_state(granularity);
		
		// insertion = max (F(i-1, j), V(i-1, j) - Wg) - Ws;
		assign_best_insertion_state(granularity);
		
		// best = max (
		// 	& dp_cell.align
		// 	& dp_cell.deletion
		// 	& dp_cell.insertion
		// );
		assign_best_best_state(granularity);
	}				
}

// 1A Upper left corner of dp matrix, should be first cell computed
void ConservidueAlignment::initialize_upper_left(AlignmentGranularity granularity) {

	// Local Smith-Waterman alignment
	if (granularity == ALIGN_LOCAL) {
		match.path_score = 0;
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score = 0;
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score = 0;
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);

	for (unsigned int i = 0; i < deletion.size(); i ++) { 
		deletion[i].path_score = 0;
		insertion[i].path_score = 0;
		deletion[i].traceback_pointer = NULL;
		insertion[i].traceback_pointer = NULL;
	}
	match.path_score = 0; // This one initialization cell does not require gaps
	
	match.traceback_pointer = NULL;
	
	if (sequence1.left_gap_factor != 0) { // Initial end gap deletion score
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			double gap_score = 
				conservidue1->gap_deletion_penalty +
				conservidue2->gap_opening_penalty +
				conservidue2->gap_open_offset(gseg); // piecewise penalty
			deletion[gseg].path_score += sequence1.left_gap_factor * gap_score;
		}
	}
	if (sequence2.left_gap_factor != 0) {
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			double gap_score = 
				conservidue2->gap_deletion_penalty +
				conservidue1->gap_opening_penalty +
				conservidue1->gap_open_offset(gseg); // piecewise penalty
			insertion[gseg].path_score += sequence2.left_gap_factor * gap_score;
		}
	}
}

// 1B Possible left end gap in sequence 1
// Top row of dp matrix
// i == 0
void ConservidueAlignment::initialize_top_row(AlignmentGranularity granularity) {
	// only need to initialize insertion path step
	// deletion score is implicitly handled in this initialization
	// there is no match to the start Conservidue
	// only insertion and match from inner matrix will query this cell
				
	// Local Smith-Waterman alignment
	if (granularity == ALIGN_LOCAL) {
		match.path_score = 0;
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score = 0;
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score = 0;
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		deletion[gseg].path_score = BAD_SCORE;
		insertion[gseg].path_score = 0;
		deletion[gseg].traceback_pointer = NULL;
	}
	match.path_score = BAD_SCORE;
	
	match.traceback_pointer = NULL;
	
	// V(0,j) = F(0,j) = -Wg1(0) - Wo2(0) - SUM[x = 1..j](We2(x))
	
	// find best predecessor of conservidue2 (only needed for non-branched alignments)
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		double best_score = BAD_SCORE;
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue2->predecessors.begin();
			 prev_res != conservidue2->predecessors.end();
			 prev_res ++) {
			int predecessor_index = prev_res->predecessor_conservidue;
			const Conservidue * previous_residue = & sequence2[predecessor_index];
			const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_residue);
		
			double test_score = 0;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->insertion[gseg].path_score;
		
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
	
		insertion[gseg].traceback_pointer = & (best_previous_cell->insertion[gseg]);
		insertion[gseg].path_score = best_score;
	
		if (sequence1.left_gap_extension_factor != 0) {
			insertion[gseg].path_score += 
			sequence1.left_gap_extension_factor * conservidue2->gap_extension_penalty(gseg);
		}
	}
}

// 1C Possible left end gap in sequence 2
// Left column of dp matrix
// j == 0
void ConservidueAlignment::initialize_left_column(AlignmentGranularity granularity) {
	// only need to initialize deletion path step
	// insertion score is implicitly handled in this initialization
	// there is no match to the start Conservidue
	// only deletion and match from inner matrix will query this cell

	// Local Smith-Waterman alignment
	if (granularity == ALIGN_LOCAL) {
		match.path_score = 0;
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score = 0;
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score = 0;
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);

	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		insertion[gseg].path_score = BAD_SCORE;
		deletion[gseg].path_score = 0;
		insertion[gseg].traceback_pointer = NULL;
	}
	match.path_score = BAD_SCORE;
	
	match.traceback_pointer = NULL;
	
	// V(i,0) = E(i,0) = -Wg2(0) - Wo1(0) - SUM[x = 1..i](We1(x))
	// TODO 
	
	// find best predecessor of conservidue1 (only needed for non-branched alignments)
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		double best_score = BAD_SCORE;
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue1->predecessors.begin();
			 prev_res != conservidue1->predecessors.end();
			 prev_res ++) {
			const Conservidue * previous_residue = & sequence1[prev_res->predecessor_conservidue];
			const ConservidueAlignment * previous_cell = get_cell(previous_residue, conservidue2);
		
			double test_score = 0;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->deletion[gseg].path_score;
		
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
	
		deletion[gseg].traceback_pointer = & (best_previous_cell->deletion[gseg]);
		deletion[gseg].path_score = best_score;
	
		if (sequence2.left_gap_extension_factor != 0) {
			deletion[gseg].path_score += 
			sequence2.left_gap_extension_factor * conservidue1->gap_extension_penalty(gseg);
		}
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
void ConservidueAlignment::assign_best_match_state(AlignmentGranularity granularity) {
	double best_score = BAD_SCORE;
	AlignmentStep * subject_step = & match;

	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	// examine all predecessor pairs
	vector<ConserviduePredecessor>::const_iterator prev_res1;
	vector<ConserviduePredecessor>::const_iterator prev_res2;
	for (prev_res1 = conservidue1->predecessors.begin();
		 prev_res1 != conservidue1->predecessors.end();
		 prev_res1 ++) {
		
		const ConserviduePredecessor & cp1 = * prev_res1;
		const Conservidue * previous_conservidue1 = & sequence1[cp1.predecessor_conservidue];
		float transition_score1 = cp1.transition_score;
		
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {

			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = & sequence2[cp2.predecessor_conservidue];
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
			for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
				test_step = & previous_cell->deletion[gseg];
				test_score = test_step->path_score + 
					transition_score1 +
					transition_score2 +
					previous_conservidue2->gap_closing_penalty * end_gap_factor2 +
					previous_conservidue2->gap_close_offset(gseg) * end_gap_factor2 +
					conservidue_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				}
			}

			// C) match to insertion
			// 	F(i-1, j-1) - Wc1(i-1) + S(i,j)
			for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
				test_step = & previous_cell->insertion[gseg];
				test_score = test_step->path_score + 
					transition_score1 +
					transition_score2 +
					previous_conservidue1->gap_closing_penalty * end_gap_factor1 +
					previous_conservidue1->gap_close_offset(gseg) * end_gap_factor1 +
					conservidue_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if match->insertion is best
			} // end for gseg on insertion
		} // end for prev_res2
	} // end for prev_res1

	// Local Smith-Waterman alignment
	if ((granularity == ALIGN_LOCAL) && (best_score <= 0)) {
		subject_step->path_score = 0;
		subject_step->traceback_pointer = NULL;
	}
}

// E(i,j) = max {
// 	E(i, j-1), 
// 	G(i, j-1) - Wg1(i) - Wo2(j-1), 
// 	F(i, j-1) - Wg1(i) - Wo2(j-1) - Wc1(i) (?)
// } - We2(j)
void ConservidueAlignment::assign_best_deletion_state(AlignmentGranularity granularity) {
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	float end_gap_factor = 1.0; // In case end gaps are free
	float end_gap_extension_factor = 1.0; // In case end gaps are free
	if (conservidue1->is_final) {
		end_gap_factor = conservidue1->parent_alignment->right_gap_factor;
		end_gap_extension_factor = conservidue1->parent_alignment->right_gap_extension_factor;
	}
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {	
		double best_score = BAD_SCORE;
		AlignmentStep * subject_step = & deletion[gseg];
	
		// Perhaps a final end gap closing penalty is needed
		float final_gap_close = 0.0;
		if (conservidue1->is_final && conservidue2->is_final) {
			final_gap_close = conservidue2->gap_closing_penalty * end_gap_factor +
				conservidue2->gap_close_offset(gseg) * end_gap_factor;
		}

		// examine all predecessors
		vector<ConserviduePredecessor>::const_iterator prev_res2;
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {
			
			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = & sequence2[cp2.predecessor_conservidue];
			float transition_score2 = cp2.transition_score;
			
			const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_conservidue2);
			double test_score;
			const AlignmentStep * test_step;
			
			// A) deletion to deletion
			// 	E(i, j-1) - We2(j)
			test_step = & previous_cell->deletion[gseg];
			test_score = test_step->path_score + 
				transition_score2 +
				final_gap_close +
				conservidue2->gap_extension_penalty(gseg) * end_gap_extension_factor;
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
				previous_conservidue2->gap_open_offset(gseg) * end_gap_factor +
				conservidue2->gap_extension_penalty(gseg) * end_gap_extension_factor;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// C) deletion to insertion
			// 	F(i, j-1) - Wg1(i) - Wo2(j-1) - Wc1(i) - We2(j)
			test_step = & previous_cell->insertion[gseg];
			test_score = test_step->path_score + 
				transition_score2 +
				final_gap_close +
				conservidue1->gap_deletion_penalty * end_gap_factor +
				previous_conservidue2->gap_opening_penalty * end_gap_factor +
				previous_conservidue2->gap_open_offset(gseg) * end_gap_factor +
				conservidue1->gap_closing_penalty + // Not part of any end gap!
				conservidue1->gap_close_offset(gseg) + // Not part of any end gap!
				conservidue2->gap_extension_penalty(gseg) * end_gap_extension_factor;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		}		
	} // end for gseg
}

// F(i,j) = max {
// 	F(i-1, j), 
// 	G(i-1, j) - Wg2(j) -Wo1(i-1), 
// 	E(i-1, j) - Wg2(j) -Wo1(i-1) - Wc2(j) (?)
// } - We1(i)
void ConservidueAlignment::assign_best_insertion_state(AlignmentGranularity granularity) {
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	
	float end_gap_factor = 1.0; // In case end gaps are free
	float end_gap_extension_factor = 1.0; // In case end gaps are free
	if (conservidue2->is_final) {
		end_gap_factor = conservidue2->parent_alignment->right_gap_factor;
		end_gap_extension_factor = conservidue2->parent_alignment->right_gap_extension_factor;
	}
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {	
		double best_score = BAD_SCORE;
		AlignmentStep * subject_step = & insertion[gseg];
	
		// Perhaps a final end gap closing penalty is needed
		float final_gap_close = 0.0;
		if (conservidue1->is_final && conservidue2->is_final) final_gap_close = conservidue1->gap_closing_penalty * end_gap_factor;
	
		// examine all predecessors
		vector<ConserviduePredecessor>::const_iterator prev_res1;
		for (prev_res1 = conservidue1->predecessors.begin();
			 prev_res1 != conservidue1->predecessors.end();
			 prev_res1 ++) {
		
			const ConserviduePredecessor & cp1 = * prev_res1;
			const Conservidue * previous_conservidue1 = & sequence1[cp1.predecessor_conservidue];
			float transition_score1 = cp1.transition_score;
		
			const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, conservidue2);
			double test_score;
			const AlignmentStep * test_step;
		
			// A) insertion to insertion
			// 	F(i-1, j) - We1(j)
			test_step = & previous_cell->insertion[gseg];
			test_score = test_step->path_score + 
				transition_score1 +
				final_gap_close +
				conservidue1->gap_extension_penalty(gseg) * end_gap_extension_factor;
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
				previous_conservidue1->gap_open_offset(gseg) * end_gap_factor +
				conservidue1->gap_extension_penalty(gseg) * end_gap_extension_factor;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// C) insertion to deletion
			// 	E(i-1, j) - Wg2(j) -Wo1(i-1) - Wc2(j) - We1(i)
			test_step = & previous_cell->deletion[gseg];
			test_score = test_step->path_score + 
				transition_score1 +
				final_gap_close +
				conservidue2->gap_deletion_penalty * end_gap_factor +
				previous_conservidue1->gap_opening_penalty * end_gap_factor +
				previous_conservidue1->gap_open_offset(gseg) * end_gap_factor +
				conservidue2->gap_closing_penalty + // This one cannot be part of end gap!
				conservidue2->gap_close_offset(gseg) + // This one cannot be part of end gap!
				conservidue1->gap_extension_penalty(gseg) * end_gap_extension_factor;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		}
	} // end for gseg
}

// best = max (
// 	& dp_cell.align
// 	& dp_cell.deletion
// 	& dp_cell.insertion
// );
void ConservidueAlignment::assign_best_best_state(AlignmentGranularity granularity) {
	best = & match;
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		if (insertion[gseg].path_score > best->path_score) best = & insertion[gseg];
		if (deletion[gseg].path_score > best->path_score) best = & deletion[gseg];
	}
}


SequenceAlignment align_profiles(
								 const SequenceAlignment & seq1, 
								 const SequenceAlignment & seq2,
								 AlignmentGranularity granularity
								 ) {
	return seq1.align(seq2, granularity);
}


