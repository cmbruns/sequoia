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

// AlignmentMethod.cpp -> ConservidueAlignment.cpp
//  routines for recursion in dynamic programming sequence alignment
//  includes class methods for classes ConservidueAlignment and AlignmentStep
//
//  $Id$
//
//  $Header$
//
//  $Log$
//  Revision 1.2  2004/06/14 16:49:57  cmbruns
//  Removed open/delete penalties from upper-left cell.  These are now handled by profile scoring Conservidue subroutines
//  Finalized switch from insert<->delete in left and top initialization.  I think Gusfield book was wrong.
//  Replace computation of end gap factors here, with conservidue based ones in new Conservidue scoring functions.
//  Made GAP_DIVERGENCE alignment behavior the default
//
//  Revision 1.1  2004/06/07 18:48:01  cmbruns
//  Divided AlignmentMethod source code files into more object releated files DPMatrix.cpp/h and ConservidueAlignment.cpp/h
//
//  Revision 1.3  2004/05/24 14:52:06  cmbruns
//  For all insertion, deletion, and divergence states, iterate over each segment of piecewise linear gap model
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

#include "ConservidueAlignment.h"
#include "DPMatrix.h"
#include "GapModel.h"
#include "Conservidue.h"
#include "SequenceAlignment.h"
// #include <vector>
// #include "AlignmentMethod.h"
// #include "AlignmentScore.h"
// #include "Exceptions.h"

AlignmentGranularity default_granularity;
double BAD_SCORE = -1e50;

AlignmentStep::AlignmentStep() : 
traceback_pointer(NULL),
path_score(MATCH_SCORE, BAD_SCORE),
parent_cell(NULL),
is_match(false)
{}

AlignmentStep::AlignmentStep(const ConservidueAlignment * parent) : 
traceback_pointer(NULL),
path_score(MATCH_SCORE, BAD_SCORE),
parent_cell(parent),
is_match(false)
{}

ostream & AlignmentStep::print_debug(ostream & os, unsigned int indent_size) const {
	string indent = "";
	for (unsigned int i=0;i<indent_size;i++) {indent += " ";}
	os << indent << "step address = " << this << endl;
	os << indent << "traceback pointer = " << traceback_pointer << endl;
	os << indent << "path_score = " << path_score << endl;
	os << indent << "parent_cell = " << parent_cell << endl;
	return os;
}


// assignment operator
ConservidueAlignment & ConservidueAlignment::operator=(const ConservidueAlignment & c2) {
	if (&c2 == this) return *this;
	
	deletion = c2.deletion;
	insertion = c2.insertion;
	match = c2.match;
	divergence = c2.divergence;
	best = c2.best;
	
	conservidue1 = c2.conservidue1;
	conservidue2 = c2.conservidue2;
	conservidue_score = c2.conservidue_score;
	
	// Correct self referential pointers
	for (unsigned int i = 0; i < deletion.size(); i ++) { 
		deletion[i].parent_cell = this;
		insertion[i].parent_cell = this;
		divergence[i].parent_cell = this;
	}
	match.parent_cell = this;
	
	for (unsigned int i = 0; i < deletion.size(); i ++) { 
		if (c2.best == &c2.insertion[i]) best = &insertion[i];
		if (c2.best == &c2.deletion[i]) best = &deletion[i];
		if (c2.best == &c2.divergence[i]) best = &divergence[i];
	}
	if (c2.best == &c2.match) best = &match;
	
	return *this;
}

ostream & ConservidueAlignment::print_debug(ostream & os, unsigned int indent_size) const {
		string indent = "";
		for (unsigned int i=0;i<indent_size;i++) {indent += " ";}
		
		os << indent << "cell address = " << this << endl;
		
		os << indent << "conservidue1 = " << conservidue1 << endl;
		os << indent << "conservidue2 = " << conservidue2 << endl;
		os << indent << "conservidue_score = " << conservidue_score << endl;
		
		os << indent << "match step = " << endl;
		match.print_debug(os, indent_size + 2);
		
		for (unsigned int i = 0; i < insertion.size(); i++) {
			os << indent << "insertion step " << i << " = " << endl;
			insertion[i].print_debug(os, indent_size + 2);
		}

		for (unsigned int i = 0; i < deletion.size(); i++) {
			os << indent << "deletion step " << i << " = " << endl;
			deletion[i].print_debug(os, indent_size + 2);
		}
		
		os << indent << "best step = " << best << endl;
		
		return os;
	}

	// default constructor
	ConservidueAlignment::ConservidueAlignment() :
	deletion(DEFAULT_GAP_SEGMENT_COUNT),
	insertion(DEFAULT_GAP_SEGMENT_COUNT),
	divergence(DEFAULT_GAP_SEGMENT_COUNT),
	match(NULL),
	best(NULL),
	conservidue1(NULL),
	conservidue2(NULL),
	conservidue_score(MATCH_SCORE, 0.0)
	{
		// Correct self referential pointers
		for (unsigned int i = 0; i < deletion.size(); i ++) { 
			deletion[i].parent_cell = this;
			insertion[i].parent_cell = this;
			divergence[i].parent_cell = this;

			deletion[i].is_match = false;
			insertion[i].is_match = false;
			divergence[i].is_match = false;
		}
		match.parent_cell = this;
		match.is_match = true;

		best = & match; // by default
	}

	// copy constructor
	ConservidueAlignment::ConservidueAlignment(const ConservidueAlignment & c2) {
		*this = c2;
	}


// ////////////////////////////// 2) Functions //////////////////////////////////////

const ConservidueAlignment * 
  ConservidueAlignment::get_cell(const Conservidue * c1, 
								 const Conservidue * c2) const {
 	int i = c1->array_sequence_index;
 	int j = c2->array_sequence_index;
 	return & (*parent_table)[i][j];
}

// One ConservidueAlignment is one cell of the dynamic programming table
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
		// TODO - move this out of recurrence and into SequenceAlignment::align
		compute_pair_score();
				
		// align = V(i-1, j-1) + S(i,j) ;
		assign_best_match_state(granularity);
		
		// deletion = max (E(i, j-1), V(i, j-1) - Wg) - Ws;
		assign_best_deletion_state(granularity);
		
		// insertion = max (F(i-1, j), V(i-1, j) - Wg) - Ws;
		assign_best_insertion_state(granularity);
		
		if (granularity.use_divergence == PAIR_DIVERGENCE)
			assign_best_divergence_state(granularity);

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
	if (granularity.align_global == false) {
		match.path_score.clear();
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score.clear();
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score.clear();
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	for (unsigned int i = 0; i < deletion.size(); i ++) { 
		deletion[i].path_score.clear();
		insertion[i].path_score.clear();
		deletion[i].traceback_pointer = NULL;
		insertion[i].traceback_pointer = NULL;
	}
	match.path_score.clear(); // This one initialization cell does not require gaps
	
	match.traceback_pointer = NULL;
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		// deletion[gseg].path_score += conservidue1->initial_new_gap_opening_score(*conservidue2, gseg);
		// insertion[gseg].path_score += conservidue2->initial_new_gap_opening_score(*conservidue1, gseg);
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
	if (granularity.align_global == false) {
		match.path_score.clear();
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score.clear();
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score.clear();
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		deletion[gseg].path_score.clear();
		deletion[gseg].traceback_pointer = NULL;
		insertion[gseg].path_score.clear();
		insertion[gseg].traceback_pointer = NULL;

		// deletion[gseg].path_score.match() = BAD_SCORE;
		insertion[gseg].path_score.match() = BAD_SCORE; // test
	}
	match.path_score.match() = BAD_SCORE;
	match.traceback_pointer = NULL;
	
	// V(0,j) = F(0,j) = -Wg1(0) - Wo2(0) - SUM[x = 1..j](We2(x))
	
	// find best predecessor of conservidue2 (only needed for non-branched alignments)
	// INSERTION
	if (0) { // test
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		AlignmentScore best_score(MATCH_SCORE, BAD_SCORE);
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue2->predecessors.begin();
			 prev_res != conservidue2->predecessors.end();
			 prev_res ++) {
			int predecessor_index = prev_res->predecessor_conservidue;
			const Conservidue * previous_residue = & sequence2[predecessor_index];
			const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_residue);
		
			AlignmentScore test_score;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->insertion[gseg].path_score;
		
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
		insertion[gseg].traceback_pointer = & (best_previous_cell->insertion[gseg]);
		insertion[gseg].path_score = best_score;
		
		insertion[gseg].path_score += conservidue2->initial_gap_extension_score(*conservidue1, gseg);

	}
	}
				
	// DELETION - test
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		AlignmentScore best_score(MATCH_SCORE, BAD_SCORE);
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue2->predecessors.begin();
			 prev_res != conservidue2->predecessors.end();
			 prev_res ++) {
			int predecessor_index = prev_res->predecessor_conservidue;
			const Conservidue * previous_residue = & sequence2[predecessor_index];
			const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_residue);
			
			AlignmentScore test_score;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->deletion[gseg].path_score;
			
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
		deletion[gseg].traceback_pointer = & (best_previous_cell->deletion[gseg]);
		deletion[gseg].path_score = best_score;
		
		deletion[gseg].path_score += conservidue1->initial_gap_extension_score(*conservidue2, gseg);
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
	if (granularity.align_global == false) {
		match.path_score.clear();
		match.traceback_pointer = NULL;
		for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
			insertion[gseg].path_score.clear();
			insertion[gseg].traceback_pointer = NULL;
			deletion[gseg].path_score.clear();
			deletion[gseg].traceback_pointer = NULL;
		}
		return;
	}
	
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);

	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		deletion[gseg].path_score.clear();
		deletion[gseg].traceback_pointer = NULL;
		insertion[gseg].path_score.clear();
		insertion[gseg].traceback_pointer = NULL;

		// insertion[gseg].path_score.match() = BAD_SCORE;
		deletion[gseg].path_score.match() = BAD_SCORE;
	}
	match.path_score.match() = BAD_SCORE;	
	match.traceback_pointer = NULL;
	
	// V(i,0) = E(i,0) = -Wg2(0) - Wo1(0) - SUM[x = 1..i](We1(x))
	
	// find best predecessor of conservidue1 (only needed for non-branched alignments)
	// DELETION
	if (0) {
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		AlignmentScore best_score(MATCH_SCORE, BAD_SCORE);
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue1->predecessors.begin();
			 prev_res != conservidue1->predecessors.end();
			 prev_res ++) {
			const Conservidue * previous_residue = & sequence1[prev_res->predecessor_conservidue];
			const ConservidueAlignment * previous_cell = get_cell(previous_residue, conservidue2);
		
			AlignmentScore test_score;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->deletion[gseg].path_score;
		
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
	
		deletion[gseg].traceback_pointer = & (best_previous_cell->deletion[gseg]);
		deletion[gseg].path_score = best_score;
	
		deletion[gseg].path_score += conservidue1->initial_gap_extension_score(*conservidue2, gseg);
	}
	}

	// INSERTION - test
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
		const ConservidueAlignment * best_previous_cell = NULL;
		AlignmentScore best_score(MATCH_SCORE, BAD_SCORE);
		vector<ConserviduePredecessor>::const_iterator prev_res;
		for (prev_res = conservidue1->predecessors.begin();
			 prev_res != conservidue1->predecessors.end();
			 prev_res ++) {
			const Conservidue * previous_residue = & sequence1[prev_res->predecessor_conservidue];
			const ConservidueAlignment * previous_cell = get_cell(previous_residue, conservidue2);
			
			AlignmentScore test_score;
			test_score += prev_res->transition_score; // zero for non-branched alignments
			test_score += previous_cell->insertion[gseg].path_score;
			
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				best_previous_cell = previous_cell;
			}
		}
		
		insertion[gseg].traceback_pointer = & (best_previous_cell->insertion[gseg]);
		insertion[gseg].path_score = best_score;
		
		insertion[gseg].path_score += conservidue2->initial_gap_extension_score(*conservidue1, gseg);
	}

}

// compute pairwise alignment score S(i,j);
void ConservidueAlignment::compute_pair_score() {
	conservidue_score.clear();

	// Contribution to score from matches
	conservidue_score += conservidue1->match_score(*conservidue2);
}

// align = V(i-1, j-1) + S(i,j) ;
// G(i,j) = max{
// 	E(i-1, j-1) - Wc2(j-1), 
// 	F(i-1, j-1) - Wc1(i-1), 
// 	G(i-1, j-1)
//  D(i-1,j-1) - (Wc2(j-1) + Wc1(i-1))/2
// } + S(i,j) 
void ConservidueAlignment::assign_best_match_state(AlignmentGranularity granularity) {
	double best_score = BAD_SCORE;
	AlignmentStep * subject_step = & match;

	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	// Perhaps a final end gap closing penalty is needed
	AlignmentScore final_gap_close;
	// TODO - merge these into a single function
	final_gap_close = 
		conservidue1->final_match_gap_closing_score(*conservidue2, 0) +
		conservidue2->final_match_gap_closing_score(*conservidue1, 0);
	AlignmentScore initial_gap_open;
	initial_gap_open = 
		conservidue1->initial_match_gap_opening_score(*conservidue2, 0) +
		conservidue2->initial_match_gap_opening_score(*conservidue1, 0);
	
	// examine all predecessor pairs
	vector<ConserviduePredecessor>::const_iterator prev_res1;
	vector<ConserviduePredecessor>::const_iterator prev_res2;
	for (prev_res1 = conservidue1->predecessors.begin();
		 prev_res1 != conservidue1->predecessors.end();
		 prev_res1 ++) {
		
		const ConserviduePredecessor & cp1 = * prev_res1;
		const Conservidue * previous_conservidue1 = & sequence1[cp1.predecessor_conservidue];
		AlignmentScore transition_score1 = cp1.transition_score;
		
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {

			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = & sequence2[cp2.predecessor_conservidue];
			AlignmentScore transition_score2 = cp2.transition_score;

			// score common to all paths
			AlignmentScore match_score = 
				transition_score1 + 
				transition_score2 + 
				conservidue_score +
				final_gap_close +
				initial_gap_open +
				// Contribution to score from gap extension
				// This is not correct for piecewise linear gap penalties, note gseg argument
				//conservidue_score += conservidue1->match_gap_extension_score(*conservidue2, 0);
				conservidue1->match_gap_extension_score(*conservidue2, 0);

			AlignmentScore match_open1_score =
				conservidue1->match_gap_opening_score(*conservidue2, 0);
			AlignmentScore match_open2_score =
				conservidue2->match_gap_opening_score(*conservidue1, 0);

			AlignmentScore match_close1_score =
				conservidue1->match_gap_closing_score(*previous_conservidue2, 0);
			AlignmentScore match_close2_score =
				conservidue2->match_gap_closing_score(*previous_conservidue1, 0);
			
			const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, previous_conservidue2);
			AlignmentScore test_score;
			const AlignmentStep * test_step;

			// A) match to match
			// G(i-1, j-1) + S(i,j) 
			test_step = & previous_cell->match;
			test_score = test_step->path_score + 
				match_score +
				match_open1_score +
				match_open2_score +
				match_close1_score +
				match_close2_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}

			// one insertion/deletion/divergence step for each segment in gap model
			for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {

				if (granularity.use_divergence == GAP_DIVERGENCE) {
					match_open1_score.clear();
					match_open2_score.clear();
				}
				
				AlignmentScore new_close_score1 = conservidue1->new_gap_closing_score(*previous_conservidue2, gseg);
				AlignmentScore new_close_score2 = conservidue2->new_gap_closing_score(*previous_conservidue1, gseg); 
				
				// B) match to deletion
				// 	E(i-1, j-1) - Wc2(j-1) + S(i,j)
				test_step = & previous_cell->deletion[gseg];
				test_score = test_step->path_score + 
					new_close_score1 +
					match_open2_score +
					match_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				}

				// C) match to insertion
				// 	F(i-1, j-1) - Wc1(i-1) + S(i,j)
				test_step = & previous_cell->insertion[gseg];
				test_score = test_step->path_score + 
					new_close_score2 +
					match_open1_score +
					match_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if match->insertion is best

				// D) match to divergence
				//  D(i-1,j-1) - (Wc2(j-1) + Wc1(i-1))/2
				if (granularity.use_divergence == PAIR_DIVERGENCE) {
					test_step = & previous_cell->divergence[gseg];
					test_score = test_step->path_score + 
						0.5 * new_close_score1 +
						0.5 * new_close_score2 +
						0.5 * match_open1_score +
						0.5 * match_open2_score +
						match_score;
					if ((test_score > best_score) || (best_score == BAD_SCORE)) {
						best_score = test_score;
						subject_step->path_score = test_score;
						subject_step->traceback_pointer = test_step;
					} // end if match->divergence is best
				} // end if pair divergence
				
			} // end for gseg
		} // end for prev_res2
	} // end for prev_res1

	// Local Smith-Waterman alignment
	if ((granularity.align_global == false) && (best_score <= 0)) {
		subject_step->path_score.clear();
		subject_step->traceback_pointer = NULL;
	}
}

// E(i,j) = max {
// 	E(i, j-1), 
// 	G(i, j-1) - Wg1(i) - Wo2(j), 
// 	F(i, j-1) - Wg1(i) - Wo2(j) - Wc1(i)
//  D(i, j-1) - Wg1(i)
// } - We2(j)
void ConservidueAlignment::assign_best_deletion_state(AlignmentGranularity granularity) {
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);
	
	// One deletion step for each segment in the piecewise linear gap model
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {	
		double best_score = BAD_SCORE;
		AlignmentStep * subject_step = & deletion[gseg];
	
		// examine all predecessors
		vector<ConserviduePredecessor>::const_iterator prev_res2;
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {
			
			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = & sequence2[cp2.predecessor_conservidue];
			AlignmentScore transition_score2 = cp2.transition_score;
			
			// A final end gap closing penalty may be needed
			AlignmentScore final_gap_close;
			final_gap_close = conservidue1->final_new_gap_closing_score(*conservidue2, gseg);
			AlignmentScore initial_gap_open;
			initial_gap_open = conservidue1->initial_new_gap_opening_score(*conservidue2, gseg);
			
			const ConservidueAlignment * previous_cell = get_cell(conservidue1, previous_conservidue2);
			AlignmentScore test_score;
			const AlignmentStep * test_step;

			// common to all scores
			AlignmentScore deletion_score = transition_score2 + 
				final_gap_close +
				initial_gap_open +
				conservidue1->new_gap_extension_score(*conservidue2, gseg);
			
			// common to gap opening scores
			AlignmentScore opening_score = deletion_score +
				conservidue1->new_gap_opening_score(*conservidue2, gseg);
			
			// A) deletion after deletion
			// 	E(i, j-1) - We2(j)
			test_step = & previous_cell->deletion[gseg];
			test_score = test_step->path_score + 
				deletion_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// B) deletion after match
			// 	G(i, j-1) - Wg1(i) - Wo2(j) - We2(j)
			test_step = & previous_cell->match;
			test_score = test_step->path_score + 
				opening_score;
			if (granularity.use_divergence != GAP_DIVERGENCE)
				test_score += conservidue2->match_gap_closing_score(*conservidue1, gseg);
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// C) deletion after insertion
			// 	F(i, j-1) - Wg1(i) - Wo2(j) - Wc1(i) - We2(j)
			// X- i
			// -X j
			test_step = & previous_cell->insertion[gseg];
			if (granularity.use_divergence == GAP_DIVERGENCE) {
				// TODO - divide up open, close, and offset scores evenly
				test_score = test_step->path_score + 
				deletion_score;
			}
			else {
				test_score = test_step->path_score + 
				conservidue2->new_gap_closing_score(*conservidue1, gseg) + 
				opening_score;
			}
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}

			// D) deletion to divergence
			//  D(i, j-1) - Wg1(i) - We2(j)
			if (granularity.use_divergence == PAIR_DIVERGENCE) {
				test_step = & previous_cell->divergence[gseg];
				test_score = test_step->path_score + 
					conservidue2->new_gap_deletion_score(*conservidue1, gseg) +
					deletion_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if deletion->divergence is best
			} // end if pair divergence

		}		
	} // end for gseg
}

// F(i,j) = max {
// 	F(i-1, j), 
// 	G(i-1, j) - Wg2(j) -Wo1(i), 
// 	E(i-1, j) - Wg2(j) -Wo1(i) - Wc2(j)
// } - We1(i)
void ConservidueAlignment::assign_best_insertion_state(AlignmentGranularity granularity) {
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	
	for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {	
		double best_score = BAD_SCORE;
		AlignmentStep * subject_step = & insertion[gseg];
	
		// Perhaps a final end gap closing penalty is needed
		AlignmentScore final_gap_close;
		final_gap_close = conservidue2->final_new_gap_closing_score(*conservidue1, gseg);
		AlignmentScore initial_gap_open;
		initial_gap_open = conservidue2->initial_new_gap_opening_score(*conservidue1, gseg);
		
		
		// examine all predecessors
		vector<ConserviduePredecessor>::const_iterator prev_res1;
		for (prev_res1 = conservidue1->predecessors.begin();
			 prev_res1 != conservidue1->predecessors.end();
			 prev_res1 ++) {
		
			const ConserviduePredecessor & cp1 = * prev_res1;
			const Conservidue * previous_conservidue1 = & sequence1[cp1.predecessor_conservidue];
			AlignmentScore transition_score1 = cp1.transition_score;
		
			const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, conservidue2);
			AlignmentScore test_score;
			const AlignmentStep * test_step;
		
			AlignmentScore insertion_score = transition_score1 + 
				final_gap_close +
				initial_gap_open + 
				conservidue2->new_gap_extension_score(*conservidue1, gseg);
			
			// common to gap opening scores
			AlignmentScore opening_score = insertion_score +
				conservidue2->new_gap_opening_score(*conservidue1, gseg);
			
			// A) insertion after insertion
			// 	F(i-1, j) - We1(j)
			test_step = & previous_cell->insertion[gseg];
			test_score = test_step->path_score + 
				insertion_score;
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// B) insertion after match
			// 	G(i-1, j) - Wg2(j) -Wo1(i-1) - We1(j)
			test_step = & previous_cell->match;
			test_score = test_step->path_score + 
				opening_score;
			if (granularity.use_divergence != GAP_DIVERGENCE)
				test_score += conservidue1->match_gap_closing_score(*conservidue2, gseg);
			if ((test_score > best_score) || (best_score == BAD_SCORE)) {
				best_score = test_score;
				subject_step->path_score = test_score;
				subject_step->traceback_pointer = test_step;
			}
		
			// C) insertion after deletion
			// We don't both insertion to deletion and deletion to insertion?
			// 	E(i-1, j) - Wg2(j) -Wo1(i-1) - Wc2(j) - We1(i)
			// Require that deletions may follow insertions, but not vice versa
			if (0) { // Don't do this
				test_step = & previous_cell->deletion[gseg];
				if (granularity.use_divergence == GAP_DIVERGENCE) {
					// TODO - divide up open, close, and offset scores evenly
					test_score = test_step->path_score + 
					insertion_score;
				}
				else {
					test_score = test_step->path_score + 
					conservidue1->new_gap_closing_score(*conservidue2, gseg) +
					opening_score;
				}
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				}
			}
			
			// D) insertion to divergence
			//  D(i-1, j) - Wg2(j) - We1(i)
			if (granularity.use_divergence == PAIR_DIVERGENCE) {
				test_step = & previous_cell->divergence[gseg];
				test_score = test_step->path_score + 
					conservidue1->new_gap_deletion_score(*conservidue2, gseg) +
					insertion_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if insertion->divergence is best
			} // end if pair divergence
			
		}
	} // end for gseg
}

// D(i,j) = max {
// 	D(i-1, j-1),
//	E(i-1, j-1), 
//	F(i-1, j-1),
//	G(i-1, j-1) - (Wo1(i) + Wo2(j))/2, 
//} - (We2(j) + We1(i)) / 2
// TODO - work out mechanics of end gaps with divergences
// The divergence state is only used in the incomplete PAIR_DIVERGENCE method
void ConservidueAlignment::assign_best_divergence_state(AlignmentGranularity granularity) {
	const SequenceAlignment & sequence1 = *(conservidue1->parent_alignment);
	const SequenceAlignment & sequence2 = *(conservidue2->parent_alignment);

	// TODO - this method is inadequate to identify all end divergences

	// examine all predecessor pairs
	vector<ConserviduePredecessor>::const_iterator prev_res1;
	vector<ConserviduePredecessor>::const_iterator prev_res2;
	for (prev_res1 = conservidue1->predecessors.begin();
		 prev_res1 != conservidue1->predecessors.end();
		 prev_res1 ++) {
		
		const ConserviduePredecessor & cp1 = * prev_res1;
		const Conservidue * previous_conservidue1 = & sequence1[cp1.predecessor_conservidue];
		AlignmentScore transition_score1 = cp1.transition_score;
		
		for (prev_res2 = conservidue2->predecessors.begin();
			 prev_res2 != conservidue2->predecessors.end();
			 prev_res2 ++) {
			
			const ConserviduePredecessor & cp2 = * prev_res2;
			const Conservidue * previous_conservidue2 = & sequence2[cp2.predecessor_conservidue];
			AlignmentScore transition_score2 = cp2.transition_score;
			
			// one insertion/deletion/divergence step for each segment in gap model
			for (unsigned int gseg = 0; gseg < deletion.size(); gseg++) {
				double best_score = BAD_SCORE;
				AlignmentStep * subject_step = & divergence[gseg];				
				
				// Perhaps a final end gap closing penalty is needed
				AlignmentScore final_gap_close;
				final_gap_close = 
					0.5 * conservidue2->final_new_gap_closing_score(*conservidue1, gseg) + 
					0.5 * conservidue1->final_new_gap_closing_score(*conservidue2, gseg);
				AlignmentScore initial_gap_open = 
					0.5 * conservidue2->initial_new_gap_opening_score(*conservidue1, gseg) + 
					0.5 * conservidue1->initial_new_gap_opening_score(*conservidue2, gseg);
				
				// score common to all paths
				AlignmentScore divergence_score = transition_score1 + 
					transition_score2 +
					final_gap_close +
					initial_gap_open +
					0.5 * conservidue1->new_gap_extension_score(*conservidue2, gseg) +
					0.5 * conservidue2->new_gap_extension_score(*conservidue1, gseg);
				
				const ConservidueAlignment * previous_cell = get_cell(previous_conservidue1, previous_conservidue2);
				AlignmentScore test_score;
				const AlignmentStep * test_step;
			
				// A) divergence to match
				//	G(i-1, j-1) - Wo1(i)/2 + Wo2(j)/2 + We2(j)/2 + We1(i)/2
				test_step = & previous_cell->match;
				test_score = test_step->path_score + 
					0.5 * conservidue2->new_gap_opening_score(*conservidue2, gseg) +
					0.5 * conservidue1->new_gap_opening_score(*conservidue2, gseg) + 
					divergence_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				}
			
				// B) divergence to deletion
				// 	E(i-1, j-1) + We2(j)/2 + We1(i)/2
				test_step = & previous_cell->deletion[gseg];
				test_score = test_step->path_score + 
					divergence_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				}
				
				// C) divergence to insertion
				// 	F(i-1, j-1) + We2(j)/2 + We1(i)/2
				test_step = & previous_cell->insertion[gseg];
				test_score = test_step->path_score + 
					divergence_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if divergence->insertion is best
				
				// D) divergence to divergence
				// 	D(i-1, j-1) + We2(j)/2 + We1(i)/2
				test_step = & previous_cell->divergence[gseg];
				test_score = test_step->path_score + 
					divergence_score;
				if ((test_score > best_score) || (best_score == BAD_SCORE)) {
					best_score = test_score;
					subject_step->path_score = test_score;
					subject_step->traceback_pointer = test_step;
				} // end if divergence->insertion is best
				
			} // end for prev_res2
		} // end for prev_res1
	} // end for gseg on insertion
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
		if ((granularity.use_divergence == PAIR_DIVERGENCE) && (divergence[gseg].path_score > best->path_score))
			best = &divergence[gseg];
	}
}


