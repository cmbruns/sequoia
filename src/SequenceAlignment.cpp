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
// Revision 1.7  2004/06/14 16:57:49  cmbruns
// Created update_begin_conservidue subroutine, but it should probably only be called in one place
// Improved logic for GAP_DIVERGENCE case in sequence_pair_score
// Created delta_sum_of_pairs_score for easier debugging in later progressive alignment steps
// Incorporated some new Conservidue members in set_weight() subroutine
//
// Revision 1.6  2004/06/04 19:30:06  cmbruns
// Updated GPL header
//
// Moved constructor methods to end of file
//
// Added initialize_from_pdb
//
// Moved all Conservidue functions to their own file Conservidue.cpp
//
// Fixed some indentation
//
// Enhanced pair score to work for local alignments, and to partially work for GAP_DIVERGENCE alignments.  TODO
//

#include "SequenceAlignment.h"
#include "MutationMatrix.h"
#include "PDBEntry.h"
// #include "ConservidueAlignment.h"
#include "DPMatrix.h"

template class map<const char, float>;

static bool dump_traceback_coordinates = false; // for debugging
static bool debug_closing = false; // dump position specific gap closing scores
static bool debug_extension = false; // dump position specific gap closing scores

// ********* SequenceAlignment methods *************

void SequenceAlignment::update_begin_conservidue() {
	if (conservidues.size() < 2) throw 99; // need at least one real conservidue
	Conservidue & begin_conservidue = conservidues[0];
	Conservidue & first_conservidue = conservidues[1];
	
	// Make most members a mirror of the first real conservidue
	begin_conservidue = first_conservidue;
	
	// Right of begin is left of first
	begin_conservidue.scaled_right_sequence_count = first_conservidue.scaled_left_sequence_count;
	begin_conservidue.scaled_extension_right_sequence_count = first_conservidue.scaled_extension_left_sequence_count;
	
	// Set left-scaled gap_deletion penalty for begin conservidue only
	if (sequence().size() == 1) {
		const BioSequence & seq = sequences[0];
		const Residue & residue = seq[0];
		begin_conservidue.scaled_gap_deletion_penalty.deletion() = residue.p_gap_deletion_penalty * weighted_sequence_count * seq.left_gap_factor;
		begin_conservidue.initial_scaled_gap_deletion_penalty.deletion() = residue.p_gap_deletion_penalty * weighted_sequence_count * seq.left_gap_factor;
		// TODO this is a kludge
		// first_conservidue.initial_scaled_gap_deletion_penalty.deletion() = 0;
	}
	
	// But maintain its special position in the sequence
	begin_conservidue.is_initial = true;
	begin_conservidue.is_final = false;
	begin_conservidue.predecessors.clear();
	begin_conservidue.array_sequence_index = 0;
	begin_conservidue.residue_counts.clear();
	begin_conservidue.residue_score_counts.assign(26, 0);
	begin_conservidue.sequence_residues.assign(sequences.size(), -1);
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
	// dp_table.print_debug(cout, 2);

	
	// ***********************************
	// **** TRACE BACK ALIGNMENT PATH ****
	// ***********************************
	
	const AlignmentStep * path_step = NULL;

	if (granularity.align_global == false) // Local Smith-Waterman alignment
		path_step = maximal_segment_pair_end;
	else if ((m > 0) && (n > 0)) path_step = dp_table[m-1][n-1].best;
	
	AlignmentScore path_score;
	if (path_step != NULL) path_score = path_step->path_score;
	// cout << "Alignment score = " << path_score << endl;
	
	// Traceback is in the opposite direction of convenient alignment construction,
	// so store coordinates for reverse iteration
	vector<long> alignment1_matches;
	vector<long> alignment2_matches;
	vector<double> debug_scores;
	AlignmentScore current_path_score = path_score;
	while (path_step != NULL) {
		int index1 = path_step->parent_cell->conservidue1->array_sequence_index;
		int index2 = path_step->parent_cell->conservidue2->array_sequence_index;
		if (dump_traceback_coordinates) cout << index1 << ", " << index2;
		if (path_step->is_match && (index1 > 0) && (index2 > 0)) {
			if (dump_traceback_coordinates) cout << "*";
			alignment1_matches.push_back(index1);
			alignment2_matches.push_back(index2);
		}
		if (dump_traceback_coordinates) {
			cout << endl;
			path_step->path_score.print_details(cout);
		}
		
		path_step = path_step->traceback_pointer;

		// Debug gap extension scores
		if (debug_extension && (path_step != NULL)) {
			AlignmentScore previous_path_score = path_step->path_score;
			double score_change = current_path_score.extension() - previous_path_score.extension();
			debug_scores.push_back(score_change);
			current_path_score = previous_path_score;
		}				
		if (debug_closing && (path_step != NULL)) {
			AlignmentScore previous_path_score = path_step->path_score;
			if (previous_path_score.closing() != current_path_score.closing()) {
				double closing_change = current_path_score.closing() - previous_path_score.closing();
				// cout << ci << " TB closing " << closing_change << endl;
			}
			current_path_score = previous_path_score;
		}				
	}
	if (debug_extension) {
		int previous_index = debug_scores.size(); // accumulate consecutive scores
		double delta_score = 0;
		for (int ci = (debug_scores.size() - 1); ci >= 0; ci --) { // each conservidue
			if ((debug_scores[ci] > 0.001) || (debug_scores[ci] < -0.001)) { // non zero extension penalty
				if (1) {															 // if ((previous_index - ci) > 1) { // accumulate consecutive scores
				// if ((previous_index - ci) > 1) { // accumulate consecutive scores
					cout << "traceback extension score ";
					cout << debug_scores.size() - previous_index + 1;
					cout << " = ";
					cout << delta_score;
					cout << endl;

					delta_score = 0;
				}
				previous_index = ci;
				delta_score += debug_scores[ci];
			}
		}
		cout << "traceback extension score ";
		cout << debug_scores.size() - previous_index;
		cout << " = ";
		cout << delta_score;
		cout << endl;
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
	
	// Now step through matches in forward direction
	answer.conservidues.clear();

	Conservidue indel_conservidue;

	int size1 = alignment1_matches.size();
	int match_index1 = alignment1_matches[size1 - 1];
	int match_index2 = alignment2_matches[size1 - 1];
	
	// Begin conservidue needs to contain open and delete penalties from the leftmost portion of the alignment
	Conservidue begin_conservidue;
	// Only if begin is a deletion will the begin_conservidue not involve seq1
	if ((alignment1_matches.size() > 0) && (match_index1 <= 1) && (match_index2 > 1)) { 
		// deletion start
		indel_conservidue = seq1[0].indel_conservidue(true);
		// indel_conservidue.initial_scaled_gap_deletion_penalty = seq1[0].initial_scaled_gap_deletion_penalty;
		begin_conservidue = indel_conservidue.combine_conservidues(seq2[0]);
	}
	else if ((alignment1_matches.size() > 0) && (match_index1 <= 1) && (match_index2 <= 1)) { 
		// match start		
		begin_conservidue = seq1[0].combine_conservidues(seq2[0]);
	}
	else { // insertion
		indel_conservidue = seq2[0].indel_conservidue(true);
		// indel_conservidue.initial_scaled_gap_deletion_penalty = seq2[0].initial_scaled_gap_deletion_penalty;
		begin_conservidue = seq1[0].combine_conservidues(indel_conservidue);
	}
	
	answer.add_conservidue(begin_conservidue);

	// Can't I make these start at 1,1, instead of at 0,0?
	unsigned long alignment1_position = 0; // Generally the next position to align in seq1
	unsigned long alignment2_position = 0; // Generally the next position to aling in seq2

	// Loop over match states, in increasing coordinate order
	for (long match = alignment1_matches.size() - 1; match >= 0; match--) {
		unsigned long match1 = alignment1_matches[match];
		unsigned long match2 = alignment2_matches[match];

		bool have_alignment1_gaps = false;
		bool have_alignment2_gaps = false;

		// Insert gaps up to match position
		if (alignment1_position < match1) {
			// Gaps will be inserted into sequence 2
			bool is_first_indel = true;
			while (alignment1_position < match1) {
				if (! seq1[alignment1_position].is_initial) {
					have_alignment2_gaps = true;
					// make sure seq2 predecessor conservidue is correct
					int pos2 = alignment2_position - 1;
					if (pos2 < 0) pos2 = 0;
					indel_conservidue = seq2[pos2].indel_conservidue(is_first_indel);
					answer.add_conservidue(seq1[alignment1_position].combine_conservidues(indel_conservidue));
					is_first_indel = false; // After first iteration, indel is no longer the first
				}
				alignment1_position ++;
			}
		}
		
		// Insert gaps up to match position
		if (alignment2_position < match2) {
			// Gaps will be inserted into sequence 1
			bool is_first_indel = true;
			while (alignment2_position < match2) {
				if (! seq2[alignment2_position].is_initial) {
					have_alignment1_gaps = true;
					// make sure seq1 predecessor conservidue is correct
					int pos1 = alignment1_position - 1;
					if (pos1 < 0) pos1 = 0;
					indel_conservidue = seq1[pos1].indel_conservidue(is_first_indel);
					Conservidue seq2_conservidue = seq2[alignment2_position];
					// seq2 gaps must precede seq1 gaps
					// if both gaps occur, adjust open close penalties for seq2 after first gaps
					if (is_first_indel && have_alignment2_gaps) {
						seq2_conservidue.scaled_gap_open_count = 0;
						seq2_conservidue.scaled_gap_deletion_penalty_opens.deletion() = 0;
						seq2_conservidue.scaled_gap_close_count = seq2_conservidue.scaled_left_sequence_count;
						seq2_conservidue.scaled_extension_gap_open_count = 0;
						seq2_conservidue.scaled_extension_gap_close_count = seq2_conservidue.scaled_extension_left_sequence_count;
					}
					answer.add_conservidue(indel_conservidue.combine_conservidues(seq2_conservidue));
					is_first_indel = false; // After first iteration, indel is no longer the first
				}
				alignment2_position ++;
			}
		}

		// insert match
		if ((! seq2[alignment2_position].is_initial) && (! seq1[alignment1_position].is_initial)) {
			Conservidue seq1_conservidue = seq1[alignment1_position];
			Conservidue seq2_conservidue = seq2[alignment2_position];
			
			// seq2 gaps must precede seq1 gaps
			// if seq1 gaps occur, adjust open close penalties for seq1 after gaps
			if (have_alignment1_gaps) {
				seq1_conservidue.scaled_gap_open_count = 0;
				seq1_conservidue.scaled_gap_deletion_penalty_opens.deletion() = 0;
				seq1_conservidue.scaled_gap_close_count = seq1_conservidue.scaled_left_sequence_count;
				seq1_conservidue.scaled_extension_gap_open_count = 0;
				seq1_conservidue.scaled_extension_gap_close_count = seq1_conservidue.scaled_extension_left_sequence_count;
			} // Only adjust seq2 if there were no gaps in seq1, otherwise it is handled in seq1 gap loop above
			else if (have_alignment2_gaps) {
				seq2_conservidue.scaled_gap_open_count = 0;
				seq2_conservidue.scaled_gap_deletion_penalty_opens.deletion() = 0;
				seq2_conservidue.scaled_gap_close_count = seq2_conservidue.scaled_left_sequence_count;
				seq2_conservidue.scaled_extension_gap_open_count = 0;
				seq2_conservidue.scaled_extension_gap_close_count = seq2_conservidue.scaled_extension_left_sequence_count;
			}

			answer.add_conservidue(seq1_conservidue.combine_conservidues(seq2_conservidue));
		}

		alignment1_position ++;
		alignment2_position ++;
	}

	bool have_alignment2_gaps = false;
	// extend gaps to end
	if (alignment1_position < m) {
		// Gaps will be inserted into sequence 2
		// Insert gaps in sequence 2
		have_alignment2_gaps = true;
		bool is_first_indel = true;
		while (alignment1_position < m) {
			if (seq1[alignment1_position].is_initial) {
				alignment1_position ++;
				continue;
			}
			int seq2_pos = alignment2_position - 1;
			if (seq2_pos < 0) seq2_pos = 0;
			indel_conservidue = seq2[seq2_pos].indel_conservidue(is_first_indel);
			
			answer.add_conservidue(seq1[alignment1_position].combine_conservidues(indel_conservidue));
			alignment1_position ++;
			is_first_indel = false;
		}
	}
	
	if (alignment2_position < n) {
		// Gaps will be inserted into sequence 1
		// Insert gaps in sequence 1
		bool is_first_indel = true;
		while (alignment2_position < n) {
			if (seq2[alignment2_position].is_initial) {
				alignment2_position ++;
				continue;
			}
			int seq1_pos = alignment1_position - 1;
			if (seq1_pos < 0) seq1_pos = 0;
			indel_conservidue = seq1[seq1_pos].indel_conservidue(is_first_indel);

			Conservidue seq2_conservidue = seq2[alignment2_position];
			// seq2 gaps must precede seq1 gaps
			// if both gaps occur, adjust open close penalties for seq2 after first gaps
			if (is_first_indel && have_alignment2_gaps) {
				seq2_conservidue.scaled_gap_open_count = 0;
				seq2_conservidue.scaled_gap_deletion_penalty_opens.deletion() = 0;
				seq2_conservidue.scaled_gap_close_count = seq2_conservidue.scaled_left_sequence_count;
				seq2_conservidue.scaled_extension_gap_open_count = 0;
				seq2_conservidue.scaled_extension_gap_close_count = seq2_conservidue.scaled_extension_left_sequence_count;
			}

			answer.add_conservidue(indel_conservidue.combine_conservidues(seq2[alignment2_position]));
			alignment2_position ++;
			is_first_indel = false;
		}
	}
	// answer.update_begin_conservidue();
		
	// answer.print_debug(); // for debugging
	
	return answer;
}

SequenceAlignment & SequenceAlignment::initialize_from_pdb_protein(const PDBChain & chain) {
	BioSequence pdb_sequence = chain.to_coord_fasta();
	initialize_from_biosequence(pdb_sequence);
	// TODO - add secondary structure and other structure information
	return *this;
}

ostream & SequenceAlignment::print_pretty(ostream & os) const {
	const SequenceAlignment & s = *this;
	// Pretty alignment output
	
	unsigned int line_length = 60;
	unsigned int conservidue_number = 0;
	unsigned int sequence_number = 0;
	
	os << "Alignment score = " << s.pair_alignment_score << endl;
	
	// Print sequence summary
	for (sequence_number = 0; sequence_number < s.sequences.size(); ++sequence_number) {
		const BioSequence & seq = s.sequences[sequence_number];
		os << sequence_number + 1 << ":\t";
		os << seq.get_id() << "\t";
		os << seq.get_title() << endl;
	}
	os << endl;
	
	// Keep track of current residue number
	vector<int> sequence_residue_numbers(s.sequences.size(), 0);
	
	// loop over line length, sequences, and conservidues
	while (conservidue_number < s.conservidues.size()) {  // number at start of block
		unsigned int block_conservidue_number = conservidue_number;
		// final fake "sequence" indicates summary row for each block
		for (sequence_number = 0; sequence_number <= s.sequences.size(); ++sequence_number) {
			
			if (sequence_number == s.sequences.size())  // Summary row
				os << "  " << " \t";
			else os << sequence_number + 1 << ":\t";
			
			const Residue * residue;
			block_conservidue_number = conservidue_number;
			for (unsigned int column = 0; column < line_length; column ++) {
				// Skip "begin" conservidues
				while ((block_conservidue_number < s.conservidues.size()) &&
					   (s.conservidues[block_conservidue_number].is_initial)) {
					block_conservidue_number ++;
				}
				if ((block_conservidue_number) >= s.conservidues.size()) break;
				const Conservidue & conservidue = s.conservidues[block_conservidue_number];
				
				// Insert blank every ten characters
				if ((column > 0) && (column % 10 == 0)) os << ' ';
				
				if (sequence_number == s.sequences.size()) { // Summary row
					if ((conservidue.residue_counts.size() == 1) &&
						(conservidue.residue_counts.begin()->second > 1.0))
						os << '*'; // identical column
					else os << ' '; // non-identical column
				}
				else { // normal sequence row, not summary
					residue = conservidue.sequence_residue(sequence_number); // doesn't work
					if (residue == NULL)  os << '-';
					else {
						os << *residue;
						if (! residue->is_gap()) // Update sequence number for end of line
							sequence_residue_numbers[sequence_number] = residue->get_residue_number();
					}
				}
				
				block_conservidue_number ++;
			}
			
			// show residue number
			if (sequence_number == s.sequences.size()) {} // Summary row
			else os << "  " << sequence_residue_numbers[sequence_number];
			
			// End of one sequence line in alignment
			os << endl;
			
		}
		
		os << endl; // break between sections
		
		conservidue_number = block_conservidue_number;
	}
	
	return os;
}

unsigned int SequenceAlignment::length() const {return conservidues.size();}
const Conservidue & SequenceAlignment::operator[](int i) const {return conservidues.at(i);}

void SequenceAlignment::set_gap_penalty(double penalty) {
	// Distribute evenly among open, close, and delete
	double each_penalty = penalty / 3.0;
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		conservidues[i].gap_opening_penalty.opening() = each_penalty;
		conservidues[i].gap_closing_penalty.closing() = each_penalty;
		conservidues[i].gap_deletion_penalty.deletion() = each_penalty;
		conservidues[i].scaled_gap_deletion_penalty.deletion() = each_penalty;
	}
}
void SequenceAlignment::set_extension_penalty(const GapModel & gap_model) {
	// Distribute evenly among open, close, and delete
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		ResidueGapParameter parameter(gap_model, conservidues[i].weighted_sequence_count);
		conservidues[i].p_gap_parameter = parameter;
	}
}
	
float SequenceAlignment::report_accuracy(const SequenceAlignment & true_alignment) {
	const SequenceAlignment & test_alignment = *this;
	// start by assuming both are pairwise alignments
		
	// index is sequence1 residue, value is sequence2 residue
	// first read matches from trusted alignment
	vector<int> true_matches(true_alignment.sequences[1].length(), -1);
	int true_match_count = 0;
	for (unsigned int i = 0; i < true_alignment.conservidues.size(); i++) {
		const Conservidue & conservidue = true_alignment.conservidues[i];
		int sequence1_residue_index = conservidue.sequence_residues[0];
		int sequence2_residue_index = conservidue.sequence_residues[1];

		if (sequence1_residue_index < 0) continue;
		if (sequence2_residue_index < 0) continue;

		const Residue & residue1 = true_alignment.sequences[0][sequence1_residue_index];
		const Residue & residue2 = true_alignment.sequences[1][sequence2_residue_index];

		if (residue1.is_gap()) continue;
		if (residue2.is_gap()) continue;

		int residue1_number = residue1.get_residue_number();
		int residue2_number = residue2.get_residue_number();

		// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
		// cout << residue1_number << ", " << residue2_number << endl;
						
		true_matches[residue1_number] = residue2_number;

		true_match_count ++;
	}
	// cout << true_match_count << " true equivalent residue pairs found" << endl;

	// next read matches from test alignment
	int true_positive_count = 0;
	int false_positive_count = 0;
	for (unsigned int i = 0; i < test_alignment.conservidues.size(); i++) {
		const Conservidue & conservidue = test_alignment.conservidues[i];
		int sequence1_residue_index = conservidue.sequence_residues[0];
		int sequence2_residue_index = conservidue.sequence_residues[1];
			
		if (sequence1_residue_index < 0) continue;
		if (sequence2_residue_index < 0) continue;
			
		const Residue & residue1 = test_alignment.sequences[0][sequence1_residue_index];
		const Residue & residue2 = test_alignment.sequences[1][sequence2_residue_index];
			
		if (residue1.is_gap()) continue;
		if (residue2.is_gap()) continue;
			
		// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
		int residue1_number = residue1.get_residue_number();
		int residue2_number = residue2.get_residue_number();
			
		if (true_matches[residue1_number] == residue2_number) {
			true_positive_count ++;
			// cout << residue1.one_letter_code() << residue2.one_letter_code() << endl;
			// cout << residue1_number << ", " << residue2_number << endl;
		}
		else {
			false_positive_count ++;
		}
		
	}
	// cout << true_positive_count << " aligned pairs correctly predicted" << endl;
	// cout << false_positive_count << " aligned pairs incorrectly predicted" << endl;

	float Q_accuracy = 100.0 * (float)true_positive_count / (float)true_match_count;
	return Q_accuracy;
}

void SequenceAlignment::add_sequence(const BioSequence & sequence) { // does not automatically update conservidues!!
	if (sequences.size() == 0)
		weighted_sequence_count = 0.0;
	sequences.push_back(sequence);
	int sequence_index = sequences.size() - 1;
	sequences.back().alignment_sequence_index = sequence_index;
	sequences.back().parent_alignment = this;
	weighted_sequence_count += sequence.get_weight();
}

void SequenceAlignment::add_sequence_automatic(const BioSequence & sequence) {
	// Unlike add_sequence, updates conservidues as well
	// Requires that the sequence length is the same as the current Alignment length
	//  (the sequence may include gaps for this)
	SequenceAlignment & alignment = *this;
	
	if (alignment.sequences.size() == 0) { // No sequences have been added yet
		SequenceAlignment new_alignment(sequence);
		alignment = new_alignment;
		// cout << new_alignment;
	}
	else {
		if (alignment.length() != (sequence.length() + 1))
			throw ALIGNMENT_LENGTH_MISMATCH_EXCEPTION();
		SequenceAlignment single_sequence_alignment(sequence);
		SequenceAlignment new_alignment(alignment);
		// 1) Add sequence structure
		new_alignment.sequences.clear();
		// From starting alignment
		for (unsigned int i = 0; i < alignment.sequences.size(); i ++)
			new_alignment.add_sequence(alignment.sequences[i]);
		// And from new sequence
		new_alignment.add_sequence(sequence);
		// 2) Add conservidues
		new_alignment.conservidues.clear();
		// TODO - make begin conservidue like the others
		Conservidue begin_conservidue;
		new_alignment.add_conservidue(begin_conservidue);  // Everyone needs a single begin
		for (unsigned int i = 1; i < alignment.length(); i ++) {
			new_alignment.add_conservidue(alignment[i].combine_conservidues(single_sequence_alignment[i]));
		}
		// new_alignment.update_begin_conservidue();
		alignment = new_alignment;
	}
}

void SequenceAlignment::add_conservidue(const Conservidue & c, unsigned int sequence_index_offset) {
	conservidues.push_back(c);
	Conservidue & conservidue = conservidues.back();
	int conservidue_index = conservidues.size() - 1;
	conservidue.array_sequence_index = conservidue_index;
	conservidue.parent_alignment = this;
	
	// TODO - I need a smarter predecessor method for POA alignments
	if (conservidue_index == 0) {
		conservidue.is_initial = true;
		
		begins.clear();
		begins.push_back(conservidue_index);
		ends.clear();
		ends.push_back(conservidue_index);
	}
	else {
		conservidue.is_initial = false;
		conservidues[conservidue_index - 1].is_final = false;
	}
	conservidue.is_final = true;
	ends[0] = conservidue_index;
	
	// Predecessor link
	ConserviduePredecessor pred;
	pred.transition_score.transition() = 0;
	pred.predecessor_conservidue = conservidue_index - 1;
	vector<ConserviduePredecessor> pred_vector(1,pred);
	conservidue.predecessors = pred_vector;
	
	// Sanity check sequence_residue mapping
	conservidue.sequence_residues.assign(sequences.size(), -1);
	for (unsigned int sr = 0; sr < c.sequence_residues.size(); sr ++) {
		if ((sr + sequence_index_offset) < conservidue.sequence_residues.size()) {
			conservidue.sequence_residues[sr + sequence_index_offset] = c.sequence_residues[sr];
		} else {
			throw 99; // something is wrong with sequence_residues logic
		}
	}
}

ostream & SequenceAlignment::print_debug(ostream & os, unsigned int indent_size) const {
	string indent = "";
	for(unsigned int i=0;i<indent_size;i++)indent += " ";
	
	os << indent << "alignment pointer = " << this << endl;
	os << indent << "left_gap_factor = " << left_gap_factor << endl;
	os << indent << "right_gap_factor = " << right_gap_factor << endl;
	os << indent << "left_gap_extension_factor = " << left_gap_extension_factor << endl;
	os << indent << "right_gap_extension_factor = " << right_gap_extension_factor << endl;
	
	os << indent << "pair_alignment_score = " << pair_alignment_score << endl;
	
	os << indent << "Conservidues:" << endl;
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		conservidues[i].print_debug(os, indent_size + 2);
		os << endl;
	}
	os << indent << "Sequences:" << endl;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		sequences[i].print_debug(os, indent_size + 2);
	}
	
	return os;
}

// Load sequence alignment from fasta file
istream & SequenceAlignment::load_fasta(istream & is) {
	SequenceAlignment & alignment = *this;
	// Read a series of sequences, and paste them together into an alignment
	while (is.good()) {
		BioSequence sequence;
		is >> sequence;
		alignment.add_sequence_automatic(sequence);
	}
	return is;
}

AlignmentScore SequenceAlignment::sequence_pair_score(unsigned int seq_index1, unsigned int seq_index2, AlignmentGranularity granularity) const {
	// Compute extension score at gap close time
	// TODO - adapt methods for different alignment granularities
	
	// Incorporate sequence weights
	double weight = sequences[seq_index1].get_weight() * sequences[seq_index2].get_weight();
	unsigned int ci;

	AlignmentScore answer;

	// Note final residue in each sequence
	const Residue * final_residue1 = NULL;
	const Residue * final_residue2 = NULL;
	for (unsigned int i = 0; i < conservidues.size(); i++) {
		const Residue * res1 = conservidues[i].sequence_residue(seq_index1);
		const Residue * res2 = conservidues[i].sequence_residue(seq_index2);
		if (res1 != NULL) final_residue1 = res1;
		if (res2 != NULL) final_residue2 = res2;
	}
	
	int gap1 = 0; // Keep track of continuty of gaps in first sequence
	int gap2 = 0;
	const Residue * latest_residue1 = NULL; // Non-gap position in sequence1
	const Residue * latest_residue2 = NULL; // Non-gap position in sequence2
	const Residue * res1 = NULL;
	const Residue * res2 = NULL;
	const BioSequence & seq1 = sequences[seq_index1];
	const BioSequence & seq2 = sequences[seq_index2];

	AlignmentScore pending_gap_penalty; // For local alignment
	bool found_match = false;

	for (ci = 0; ci < conservidues.size(); ci++) {
		res1 = conservidues[ci].sequence_residue(seq_index1);
		res2 = conservidues[ci].sequence_residue(seq_index2);
		
		if ((res1 == NULL) && (res2 == NULL)) continue; // Both are gaps, no effect on score

		// MATCH
		if ((res1 != NULL) && (res2 != NULL)) { // if match
			answer.match() += blosum62.get_score(res1->one_letter_code(), res2->one_letter_code());

			double end_gap_factor1 = 1.0;
			double end_gap_extension_factor1 = 1.0;
			if (latest_residue1 == NULL) { // closing left gap in sequence1
				end_gap_factor1 = seq1.left_gap_factor;
				end_gap_extension_factor1 = seq1.left_gap_extension_factor;
			}
			double end_gap_factor2 = 1.0;
			double end_gap_extension_factor2 = 1.0;
			if (latest_residue2 == NULL) { // closing left gap in sequence2
				end_gap_factor2 = seq2.left_gap_factor;
				end_gap_extension_factor2 = seq2.left_gap_extension_factor;
			}
			
			// Gap closing penalty
			if (granularity.use_divergence == GAP_DIVERGENCE) {
				// TODO - this is not quite correct, using seq2 parameters only
				// Seq2 to is the one with the correct begin gap factors though
				// TODO - end gap factors must be made more consistent
				double delta_extension = seq2.p_gap_model->total_penalty(gap1+gap2) * end_gap_extension_factor2 * end_gap_extension_factor1;
				pending_gap_penalty.extension() += delta_extension;
				if (debug_extension && delta_extension != 0) 
					cout << ci << " extension " << delta_extension << endl;
			}
			if (gap1 > 0) {
				// Only one of the two gaps need this condition
				if ((granularity.use_divergence != GAP_DIVERGENCE) || (gap2 == 0)) {
					double closing_change = res2->p_gap_closing_penalty * end_gap_factor1;
					if (debug_closing) // debug closing
						if (closing_change != 0) cout << ci << " closing " << closing_change << endl;
					pending_gap_penalty.closing() += closing_change;
				}
				if (granularity.use_divergence != GAP_DIVERGENCE)	{	
					double delta_extension = seq2.p_gap_model->total_penalty(gap1) * end_gap_extension_factor1;
					pending_gap_penalty.extension() += delta_extension;
					if (debug_extension && delta_extension != 0) 
						cout << ci << " extension " << delta_extension << endl;
				}
			}
			if (gap2 > 0) {
				double closing_change = res1->p_gap_closing_penalty * end_gap_factor2;
				if (debug_closing) // debug closing
					if (closing_change != 0) cout << ci << " closing " << closing_change << endl;
				pending_gap_penalty.closing() += closing_change;
				if (granularity.use_divergence != GAP_DIVERGENCE) {
					double delta_extension = seq1.p_gap_model->total_penalty(gap2) * end_gap_extension_factor2;
					pending_gap_penalty.extension() += delta_extension;
					if (debug_extension && delta_extension != 0) 
						cout << ci << " extension " << delta_extension << endl;
				}
			}
			
			gap1 = 0;
			gap2 = 0;

			// Be careful about end gap penalties in local alignment
			if ((!found_match) && (!granularity.align_global)) {} // Don't count initial end gaps in local alignment
			else // This is not the first match in a local alignment
				answer += pending_gap_penalty; // They are not final end gaps
			pending_gap_penalty.clear();

			found_match = true;			
		} // end if match

		// Insert
		else if (res1 == NULL) { // Gap in sequence 1

			double end_gap_factor = 1.0;
			if (latest_residue1 == NULL) end_gap_factor = seq1.left_gap_factor;
			
			if (latest_residue1 == final_residue1) // final residue end gap
				end_gap_factor = seq1.right_gap_factor;
			
			// gap opening penalty
			if ((gap1 == 0) &&
				((granularity.use_divergence != GAP_DIVERGENCE) || 
				 (gap2 == 0) // don't count both opens with gap divergence model
				)
			   )
			{
				pending_gap_penalty.opening() += res2->p_gap_opening_penalty * end_gap_factor;
				if (latest_residue1 != NULL) // not Initial end gap
					pending_gap_penalty.deletion() += latest_residue1->p_gap_deletion_penalty * end_gap_factor; // TODO final end factor
				else // Initial end gap
					pending_gap_penalty.deletion() += DEFAULT_GAP_DELETION_PENALTY * end_gap_factor; // TODO final end factor
			}
			
			gap1 ++;
			// gap2 = 0;
		}
		
		// Delete
		else if (res2 == NULL) { // Gap in sequence 2

			double end_gap_factor = 1.0;
			if (latest_residue2 == NULL) end_gap_factor = seq2.left_gap_factor;
			
			if (latest_residue2 == final_residue2) // final residue end gap
				end_gap_factor = seq2.right_gap_factor;

			// gap opening penalty
			if ((gap2 == 0) &&
				((granularity.use_divergence != GAP_DIVERGENCE) || 
				 (gap1 == 0) // don't count both opens with gap divergence model
				 )
				)
			{
				pending_gap_penalty.opening() += res1->p_gap_opening_penalty * end_gap_factor;
				if (latest_residue2 != NULL) // not Initial end gap
					pending_gap_penalty.deletion() += latest_residue2->p_gap_deletion_penalty * end_gap_factor;
				else // Initial end gap
					pending_gap_penalty.deletion() += DEFAULT_GAP_DELETION_PENALTY * seq2.left_gap_factor;
			}
			
			gap2 ++;
			// gap1 = 0;
		}

		if (res1 != NULL) latest_residue1 = res1;
		if (res2 != NULL) latest_residue2 = res2;
		
	} // end for conservidue
	
	// closing end gaps
	if (granularity.use_divergence == GAP_DIVERGENCE) {
		// TODO - this is not quite correct, using seq1 parameters only
		double delta_extension = seq1.p_gap_model->total_penalty(gap1+gap2) * seq1.right_gap_extension_factor;
		pending_gap_penalty.extension() += delta_extension;
		if (debug_extension && delta_extension != 0) 
			cout << ci << " extension " << delta_extension << endl;
	}
	if ((gap1 > 0) && (latest_residue2 != NULL)) 
	{ // Terminal gap in seq1
		// Only one of the two gaps need this condition
		if ((granularity.use_divergence != GAP_DIVERGENCE) || (gap2 == 0)) {
			double closing_change = latest_residue2->p_gap_closing_penalty * seq1.right_gap_factor;
			if (debug_closing) // debug closing
				if (closing_change != 0) cout << ci << " closing " << closing_change << endl;
			pending_gap_penalty.closing() += closing_change;
		}
		if (granularity.use_divergence != GAP_DIVERGENCE) {
			double delta_extension = seq2.p_gap_model->total_penalty(gap1) * seq1.right_gap_extension_factor;
			pending_gap_penalty.extension() += delta_extension;
			if (debug_extension && delta_extension != 0) 
				cout << ci << " extension " << delta_extension << endl;
		}
	}
	if ((gap2 > 0) && (latest_residue1 != NULL)) 
	{ // Terminal gap in seq2
		double closing_change = latest_residue1->p_gap_closing_penalty * seq2.right_gap_factor;
		if (debug_closing) // debug closing
			if (closing_change != 0) cout << ci << " closing " << closing_change << endl;
		pending_gap_penalty.closing() += closing_change;
		if (granularity.use_divergence != GAP_DIVERGENCE) {
			double delta_extension = seq1.p_gap_model->total_penalty(gap2) * seq2.right_gap_extension_factor;
			pending_gap_penalty.extension() += delta_extension;
			if (debug_extension && delta_extension != 0) 
				cout << ci << " extension " << delta_extension << endl;
		}
	}
	
	if (granularity.align_global) answer += pending_gap_penalty;
	pending_gap_penalty.clear();
	
	answer *= weight;
	
	return answer;
}

AlignmentScore SequenceAlignment::sum_of_pairs_score(AlignmentGranularity granularity) const {
	AlignmentScore answer;
	for (unsigned int i = 0; i < sequences.size(); i++) {
		for (unsigned int j = i + 1; j < sequences.size(); j++) {
			AlignmentScore pair_score = sequence_pair_score(i,j,granularity);
			answer += pair_score;
			// cout << "Pair score (" << i + 1 << ", " << j + 1 << ") = " << pair_score << endl;
		}
	}
	return answer;
}

// Compute sum of pairs score corresponding to latest alignment step
// border_sequence is the index of the first sequence in the second half of the alignment
AlignmentScore SequenceAlignment::delta_sum_of_pairs_score(unsigned int border_sequence, AlignmentGranularity granularity) const {
	AlignmentScore answer;
	for (unsigned int i = 0; i < border_sequence; i++) {
		for (unsigned int j = border_sequence; j < sequences.size(); j++) {
			AlignmentScore pair_score = sequence_pair_score(i,j,granularity);
			answer += pair_score;
			// cout << "Pair score (" << i + 1 << ", " << j + 1 << ") = " << pair_score << endl;
		}
	}
	return answer;
}

SequenceAlignment & SequenceAlignment::initialize_from_biosequence(const BioSequence & seq0) {
	sequences.clear();
	add_sequence(seq0);
	BioSequence & seq = sequences.back(); // Use actual stored sequence copy from now on
	int sequence_index = sequences.size() - 1;
	weighted_sequence_count = seq0.get_weight();
	
	conservidues.clear(); // delete previous contents
	begins.clear();
	ends.clear();
	
	// Add "begin" conservidue - does not actually contain residues (but end does!)
	Conservidue begin_conservidue;
	add_conservidue(begin_conservidue);
	Conservidue * current_conservidue_pointer = & conservidues.back();
	
	// Loop over residues in the sequence
	for (unsigned int i = 0; i < seq.length(); ++i) {
		const Residue & residue = seq[i];
		Conservidue current_conservidue0(residue); // initialize from single residue
		
		add_conservidue(current_conservidue0);
		Conservidue & current_conservidue = conservidues.back();
		
		if (i == 0) {
			update_begin_conservidue();
		}
		
		// store sequence/residue relationship in conservidue
		current_conservidue.sequence_residues[sequence_index] = i;
		
		// Set up links for linear sequence
		ConserviduePredecessor back_link;
		back_link.predecessor_conservidue = i;
		back_link.transition_score.transition() = 0; // no penalty for following normal sequence
		current_conservidue.predecessors.clear();
		current_conservidue.predecessors.push_back(back_link); // Link to previous conservidue
		
		// Update for next round
		current_conservidue_pointer = & current_conservidue; // Will this pointer remain valid?		
	}
	
	return *this;
}

SequenceAlignment & SequenceAlignment::operator=(const SequenceAlignment & alignment2) {
	if (this == &alignment2) return *this;
	begins = alignment2.begins;
	ends = alignment2.ends;
	left_gap_factor = alignment2.left_gap_factor;
	right_gap_factor = alignment2.right_gap_factor;
	left_gap_extension_factor = alignment2.left_gap_extension_factor;
	right_gap_extension_factor = alignment2.right_gap_extension_factor;
	pair_alignment_score = alignment2.pair_alignment_score;
	weighted_sequence_count = alignment2.weighted_sequence_count;
	
	sequences.clear();
	for (unsigned int i = 0; i < alignment2.sequences.size(); i++) {
		add_sequence(alignment2.sequences[i]);
	}
	
	conservidues.clear();
	for (unsigned int i = 0; i < alignment2.conservidues.size(); i++) {
		add_conservidue(alignment2.conservidues[i]);
	}
	
	return *this;
}

SequenceAlignment::SequenceAlignment() : 
left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
weighted_sequence_count(0), 
pair_alignment_score()
{} // simple constructor
   // need assignment operator and copy constructor so that conservidue.parent_alignment pointers are properly set

SequenceAlignment::SequenceAlignment(const SequenceAlignment & alignment2) {
	*this = alignment2;
}

SequenceAlignment::SequenceAlignment(const BioSequence & seq0) :
left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
pair_alignment_score () 
{
	initialize_from_biosequence(seq0);
}

SequenceAlignment::SequenceAlignment(const char * seq_string)  :
left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
pair_alignment_score () 
{
	BioSequence biosequence(seq_string);
	initialize_from_biosequence(biosequence);
}

// Wrapper for alignment member function
SequenceAlignment align_profiles(
								 const SequenceAlignment & seq1, 
								 const SequenceAlignment & seq2,
								 AlignmentGranularity granularity
								 ) {
	return seq1.align(seq2, granularity);
}


istream & operator>>(istream & is, SequenceAlignment & alignment) {
	alignment.load_fasta(is);
	return is;
}

// simple pretty alignment output routine
ostream & operator<<(ostream & os, const SequenceAlignment & s) {
	s.print_pretty(os);
	return os;
}



