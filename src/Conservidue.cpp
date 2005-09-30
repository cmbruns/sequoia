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
// Revision 1.3  2004/06/14 16:42:19  cmbruns
// Numerous new members and functions to handle details of profile gap penalties, toward matching sum of pairs score
//
// Revision 1.2  2004/06/04 19:04:40  cmbruns
// Initial version of indel_conservidue subroutine, intended to help clean up the trace back code in AlignmentMethod.cpp
//
// Implemented residue_score_count member, to precompute scores for faster alignment
//
// Revision 1.1  2004/05/29 22:40:21  cmbruns
// Created separate files for Conservidue, separate from SequenceAlignment
//

#include "Conservidue.h"
#include "SequenceAlignment.h"
#include "MutationMatrix.h"

// A single residue or a set of aligned residues in a sequence alignment.
// Corresponds to one column of a multiple sequence alignment

// ********* Conservidue methods *************

// Precomputed score vector for efficient match score computation
float Conservidue::get_residue_score_count(char amino_acid) const {
	int index = toupper(amino_acid - 'A');
	if (index < 0) return 0;
	if (index > ('Z' - 'A')) return 0;
	return residue_score_counts[index];
}
float & Conservidue::residue_score_count(char amino_acid) {
	int index = toupper(amino_acid - 'A');
	if (index < 0) throw 99;
	if (index > ('Z' - 'A')) throw 99;
	return residue_score_counts[index];
}

// is_first_indel communicates whether this is the first indel column after the subject conservidue
Conservidue Conservidue::indel_conservidue(bool is_first_indel) const {
	Conservidue answer;
	
	const Conservidue * before = this;
	const Conservidue * after = NULL;
	if (!before->is_final)
		after = &(before->parent_alignment->operator[](before->array_sequence_index + 1));
	
	if ((before->is_initial) && (after == NULL)) throw 99;
	
	const SequenceAlignment * sequence_alignment = NULL;
	
	if (before->is_initial) { // Left end gap
		sequence_alignment = after->parent_alignment;
	}
	else if (after == NULL) { // Right end gap
		sequence_alignment = before->parent_alignment;

		// new way
		answer.scaled_extension_gap_count = before->scaled_extension_gap_count +
			before->scaled_extension_right_sequence_count;
		answer.scaled_gap_count = before->scaled_gap_count +
			before->scaled_right_sequence_count;
		// This is the one place to populate final_scaled_gap_count with non-zero
		answer.final_scaled_gap_count = answer.scaled_gap_count;
		answer.final_scaled_extension_gap_count = answer.scaled_extension_gap_count;
	}
	else { // indel between two matches
		sequence_alignment = before->parent_alignment;

		// new way
		answer.scaled_extension_gap_count = before->scaled_extension_gap_count +
			before->scaled_extension_right_sequence_count;
		answer.scaled_gap_count = before->scaled_gap_count +
			before->scaled_right_sequence_count;		
	}
	
	// There are no sequence residues in an indel conservidue
	answer.weighted_sequence_count = 0;
	answer.weighted_gap_count = sequence_alignment->weighted_sequence_count;

	if (before->is_initial) { // Left end gap
		answer.scaled_extension_gap_count = after->scaled_extension_gap_count +
		after->scaled_extension_left_sequence_count;
		answer.scaled_gap_count = after->scaled_gap_count +
			after->scaled_left_sequence_count;
		// This is the one place to populate initial_scaled_gap_count with non-zero
		answer.initial_scaled_gap_count = answer.scaled_gap_count;
		answer.initial_scaled_extension_gap_count = answer.scaled_extension_gap_count;
		answer.initial_gap_count = answer.weighted_gap_count;
	}
	
	// gap open count depends upon whether before conservidue is EXACTLY before
	answer.scaled_gap_open_count = 0;
	answer.scaled_extension_gap_open_count = 0;
	answer.scaled_gap_deletion_penalty.clear();
	answer.scaled_gap_deletion_penalty_opens.clear();
	// answer.initial_scaled_gap_deletion_penalty.clear();
	// answer.initial_scaled_gap_deletion_penalty_opens.clear();
	if (is_first_indel) {
		answer.scaled_gap_open_count = before->scaled_right_sequence_count;
		answer.scaled_extension_gap_open_count = before->scaled_extension_right_sequence_count;
		answer.scaled_gap_deletion_penalty_opens = before->scaled_gap_deletion_penalty;
		// answer.initial_scaled_gap_deletion_penalty_opens = before->initial_scaled_gap_deletion_penalty;
	}
	
	// Initial end gaps are treated as gap opens
	if (is_first_indel && before->is_initial) {
		answer.scaled_gap_open_count = before->scaled_left_sequence_count;
		answer.scaled_extension_gap_open_count = before->scaled_extension_left_sequence_count;
		// the begin conservidue has left scaled scaled_gap_deletion_penalty
		// answer.scaled_gap_deletion_penalty_opens = before->scaled_gap_deletion_penalty;  // redundant
		// answer.initial_scaled_gap_deletion_penalty = after->initial_scaled_gap_deletion_penalty;
	}
	answer.initial_scaled_gap_deletion_penalty = before->initial_scaled_gap_deletion_penalty;
	
	// Pure indels can not be gap close positions
	// (the position following gap can be gap close)
	answer.scaled_gap_close_count = 0;
	answer.scaled_extension_gap_close_count = 0;
	
	answer.scaled_extension_left_sequence_count = 0;
	answer.scaled_left_sequence_count = 0;
	answer.scaled_extension_right_sequence_count = 0;
	answer.scaled_right_sequence_count = 0;
	answer.final_gap_closing_penalty.closing() = 0;
	answer.initial_gap_opening_penalty.opening() = 0;
	
	answer.sequence_residues.assign(sequence_alignment->sequences.size(), -1);
	
	return answer;
}

// Simple sum of pairs formula
// A(x,y) = match score = SUM(over i,j) [res(i,x) * res(j,y) * Score(i,j)]
// Requires time proportional to the number of residue types in the profile with
// the fewest residue types
AlignmentScore Conservidue::match_score(const Conservidue & c2) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;
	
	// For maximum efficieny, loop over the conservidue with fewer residue types
	const Conservidue * count_conservidue = NULL;
	const Conservidue * score_conservidue = NULL;
	if (c1.residue_counts.size() > c2.residue_counts.size()) { // Use smaller c2
		count_conservidue = &c2;
		score_conservidue = &c1;
	}
	else { // Use smaller c1
		count_conservidue = &c1;
		score_conservidue = &c2;
	}

	// accumulate pair score using residue_score_counts method
	map<char,float>::const_iterator counts;
	for (counts = count_conservidue->residue_counts.begin(); 
		 counts != count_conservidue->residue_counts.end();
		 counts ++) {
		char count_residue = counts->first;
		float count_weight = counts->second;
		answer.match() += count_weight * score_conservidue->get_residue_score_count(count_residue);
	}

	return answer;
}

AlignmentScore Conservidue::initial_gap_extension_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;
	
	// use precomputed factors
	answer = 		
		(c1.scaled_extension_gap_count + c1.scaled_extension_left_sequence_count) * c2.gap_extension_penalty(gseg) +
		c1.initial_new_gap_opening_score(c2, gseg);
	
	return answer;
}

// C(x,y) = match gap extension score = EP(y) * gap_count(x) + EP(x) * gap_count(y)
AlignmentScore Conservidue::match_gap_extension_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;

	// use precomputed factors
	answer = c1.gap_extension_penalty(gseg) * 
		c2.scaled_extension_gap_count +
		c2.gap_extension_penalty(gseg) *
		c1.scaled_extension_gap_count;
	
	return answer;
}

// B(x,y) = new gap extension score = EP(y) * right_gap_count(x)?
AlignmentScore Conservidue::new_gap_extension_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;
	
	// use precomputed factors
	answer = 
		c2.gap_extension_penalty(gseg) *
		(c1.scaled_extension_gap_count + c1.scaled_extension_right_sequence_count);
	
	return answer;
}

AlignmentScore Conservidue::initial_gap_opening_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;

	answer = c2.gap_opening_penalty * c1.scaled_left_sequence_count +
		c2.gap_open_offset(gseg) * c1.scaled_extension_left_sequence_count +
		c1.scaled_gap_deletion_penalty * c2.weighted_sequence_count;
				
	return answer;
}

// D(x,y) = match gap open score = open_count(x) * OP(y)
// H(x,y) = match gap deletion score = open_count(x) * DP * seq_count(y) ???? (no pos specific?, or store open_DP(x) as well...)
AlignmentScore Conservidue::match_gap_opening_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;
	
	answer = c2.gap_opening_penalty * c1.scaled_gap_open_count +
		c2.gap_open_offset(gseg) * c1.scaled_extension_gap_open_count +
		c1.scaled_gap_deletion_penalty_opens * c2.weighted_sequence_count;
	
	return answer;
}

// E(x,y) = new gap open score = seq_count(x) * OP(y)
// J(x,y) = new gap deletion score = seq_count(x) * DP * seq_count(y) = DP(x) * seq_count(y)
AlignmentScore Conservidue::new_gap_opening_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;

	answer = c1.scaled_right_sequence_count * c2.gap_opening_penalty +
		c1.scaled_extension_right_sequence_count * c2.gap_open_offset(gseg) +
		c2.new_gap_deletion_score(c1, gseg);
				
	return answer;
}

// J(x,y) = new gap deletion score = seq_count(x) * DP * seq_count(y) = DP(x) * seq_count(y)
AlignmentScore Conservidue::new_gap_deletion_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;
	
	answer = c1.weighted_sequence_count * c2.scaled_gap_deletion_penalty;
				
	return answer;
}

// F(x,y) = match gap close score = close_count(x) * CP(y)
AlignmentScore Conservidue::match_gap_closing_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;
	AlignmentScore answer;
	
	// debug
	// double score1 = c2.gap_closing_penalty;
	// double count1 = c1.scaled_gap_close_count;
	// double score2 = c2.gap_close_offset(gseg);
	// double count2 = c1.scaled_extension_gap_close_count;
	
	answer = c2.gap_closing_penalty * c1.scaled_gap_close_count +
		c2.gap_close_offset(gseg) * c1.scaled_extension_gap_close_count;
	
	return answer;
}

// G(x,y) = new gap close score = seq_count(x) * CP(y)
AlignmentScore Conservidue::new_gap_closing_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;

	answer = c2.gap_closing_penalty * c1.scaled_left_sequence_count +
		c2.gap_close_offset(gseg) * c1.scaled_extension_left_sequence_count;
				
	return answer;
}


// Special case of closing a final end gap
AlignmentScore Conservidue::final_match_gap_closing_score(const Conservidue & c2, unsigned int gseg) const {
	return this->final_scaled_gap_count * c2.final_gap_closing_penalty +
	// TODO - the penalty part of this terminal extension offset is still to big
		   this->final_scaled_extension_gap_count * c2.gap_close_offset(gseg);
}
AlignmentScore Conservidue::final_new_gap_closing_score(const Conservidue & c2, unsigned int gseg) const {
	return (this->final_scaled_sequence_count + this->final_scaled_gap_count) * c2.final_gap_closing_penalty + 
		   (this->final_scaled_extension_sequence_count + this->final_scaled_extension_gap_count) * c2.gap_close_offset(gseg);
}
AlignmentScore Conservidue::initial_match_gap_opening_score(const Conservidue & c2, unsigned int gseg) const {
	// debug
	double count1 = this->initial_scaled_gap_count;
	double score1 = c2.initial_gap_opening_penalty;
	double count3 = c2.initial_sequence_count;
	double score3 = this->initial_scaled_gap_deletion_penalty;

	return 
	this->initial_scaled_gap_count * c2.initial_gap_opening_penalty +
	this->initial_scaled_extension_gap_count * c2.gap_open_offset(gseg) +	
	// TODO - include deletion score component
	c2.initial_sequence_count * this->initial_scaled_gap_deletion_penalty;
}
AlignmentScore Conservidue::initial_new_gap_opening_score(const Conservidue & c2, unsigned int gseg) const {
	// debug
	double count1 = this->initial_scaled_sequence_count + this->initial_scaled_gap_count;
	double score1 = c2.initial_gap_opening_penalty;
	double count3 = c2.initial_sequence_count;
	double score3 = this->initial_scaled_gap_deletion_penalty;

	return 
	(this->initial_scaled_sequence_count + this->initial_scaled_gap_count) * c2.initial_gap_opening_penalty + 
	(this->initial_scaled_extension_sequence_count + this->initial_scaled_extension_gap_count) * c2.gap_open_offset(gseg) +
	// TODO - include deletion score component
	c2.initial_sequence_count * this->initial_scaled_gap_deletion_penalty;
}


const Residue * Conservidue::sequence_residue(unsigned int sequence_index) const {
	int residue_index = sequence_residues[sequence_index];
	if (residue_index < 0) return NULL; // sequence not found
	const BioSequence & seq = parent_alignment->sequences[sequence_index];
	const Residue * res_ptr = & seq[residue_index];
	return res_ptr;
}

// Create a combined conservidue from the alignment of two matching conservidues
Conservidue Conservidue::combine_conservidues(const Conservidue & conservidue2) const {
	const Conservidue & conservidue1 = *this;
	Conservidue answer = conservidue1;
	
	answer.weighted_sequence_count += conservidue2.weighted_sequence_count;
	answer.weighted_gap_count += conservidue2.weighted_gap_count;

	answer.gap_opening_penalty += conservidue2.gap_opening_penalty;
	answer.gap_closing_penalty += conservidue2.gap_closing_penalty;
	answer.p_gap_parameter = conservidue1.p_gap_parameter.combine(conservidue2.p_gap_parameter);
	answer.gap_deletion_penalty += conservidue2.gap_deletion_penalty;
	answer.parent_alignment = NULL;
	
	answer.scaled_extension_gap_count += conservidue2.scaled_extension_gap_count;
	answer.scaled_gap_count += conservidue2.scaled_gap_count;
	answer.scaled_gap_open_count += conservidue2.scaled_gap_open_count;
	answer.scaled_gap_close_count += conservidue2.scaled_gap_close_count;
	answer.scaled_extension_left_sequence_count += conservidue2.scaled_extension_left_sequence_count;
	answer.scaled_left_sequence_count += conservidue2.scaled_left_sequence_count;
	answer.scaled_extension_right_sequence_count += conservidue2.scaled_extension_right_sequence_count;
	answer.scaled_right_sequence_count += conservidue2.scaled_right_sequence_count;

	answer.scaled_extension_gap_open_count += conservidue2.scaled_extension_gap_open_count;
	answer.scaled_extension_gap_close_count += conservidue2.scaled_extension_gap_close_count;
	answer.final_scaled_gap_count += conservidue2.final_scaled_gap_count;
	answer.final_scaled_sequence_count += conservidue2.final_scaled_sequence_count;
	answer.final_scaled_extension_gap_count += conservidue2.final_scaled_extension_gap_count;
	answer.final_scaled_extension_sequence_count += conservidue2.final_scaled_extension_sequence_count;
	answer.initial_scaled_gap_count += conservidue2.initial_scaled_gap_count;
	answer.initial_scaled_sequence_count += conservidue2.initial_scaled_sequence_count;
	answer.initial_scaled_extension_gap_count += conservidue2.initial_scaled_extension_gap_count;
	answer.initial_scaled_extension_sequence_count += conservidue2.initial_scaled_extension_sequence_count;
	answer.scaled_gap_deletion_penalty += conservidue2.scaled_gap_deletion_penalty;
	answer.scaled_gap_deletion_penalty_opens += conservidue2.scaled_gap_deletion_penalty_opens;
	answer.final_gap_closing_penalty += conservidue2.final_gap_closing_penalty;
	answer.initial_gap_opening_penalty += conservidue2.initial_gap_opening_penalty;
	
	answer.initial_sequence_count += conservidue2.initial_sequence_count;
	answer.initial_gap_count += conservidue2.initial_gap_count;
	answer.initial_scaled_gap_deletion_penalty += conservidue2.initial_scaled_gap_deletion_penalty;
	// answer.initial_scaled_gap_deletion_penalty_opens += conservidue2.initial_scaled_gap_deletion_penalty_opens;

	// Fold in sequence to residues hash
	for (unsigned int i = 0; i < conservidue2.sequence_residues.size(); i++) {
		answer.sequence_residues.push_back(conservidue2.sequence_residues[i]);
	}
	
	// Fold in residue type counts hash
	map<char, float>::const_iterator res_count2;
	for (res_count2 = conservidue2.residue_counts.begin();
		 res_count2 != conservidue2.residue_counts.end();
		 res_count2 ++) {
		answer.residue_counts[res_count2->first] += res_count2->second;
	}

	// Fold in residue score counts hash
	vector<float>::const_iterator res_score_count2;
	for (unsigned int r = 0; r < residue_score_counts.size(); r++) {
		answer.residue_score_counts[r] += conservidue2.residue_score_counts[r];
	}
	
	return answer;
}

AlignmentScore Conservidue::gap_extension_penalty(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.extension_penalty(gseg);
}
AlignmentScore Conservidue::gap_open_offset(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.open_offset(gseg);
}
AlignmentScore Conservidue::gap_close_offset(unsigned int gseg) const { // -log2 probability of this residue adding to a loop
	return p_gap_parameter.close_offset(gseg);
}

ostream & Conservidue::print_debug(ostream & os, unsigned int indent_size) const {
	string indent = "";
	for(unsigned int i=0;i<indent_size;i++)indent += " ";
	
	os << indent << "conservidue pointer = " << this << endl;
	os << indent << "parent_alignment = " << parent_alignment << endl;
	os << indent << "array_sequence_index = " << array_sequence_index << endl;
	
	os << indent << "weighted_sequence_count = " << weighted_sequence_count << endl;
	os << indent << "weighted_gap_count = " << weighted_gap_count << endl;

	os << indent << "scaled_extension_gap_count = " << scaled_extension_gap_count << endl;
	os << indent << "scaled_gap_count = " << scaled_gap_count << endl;
	os << indent << "scaled_gap_open_count = " << scaled_gap_open_count << endl;
	os << indent << "scaled_gap_close_count = " << scaled_gap_close_count << endl;
	os << indent << "scaled_extension_left_sequence_count = " << scaled_extension_left_sequence_count << endl;
	os << indent << "scaled_left_sequence_count = " << scaled_left_sequence_count << endl;
	os << indent << "scaled_extension_right_sequence_count = " << scaled_extension_right_sequence_count << endl;
	os << indent << "scaled_right_sequence_count = " << scaled_right_sequence_count << endl;

	os << indent << "scaled_extension_gap_open_count = " << scaled_extension_gap_open_count << endl;
	os << indent << "scaled_extension_gap_close_count = " << scaled_extension_gap_close_count << endl;
	os << indent << "final_scaled_gap_count = " << final_scaled_gap_count << endl;
	os << indent << "final_scaled_sequence_count = " << final_scaled_sequence_count << endl;
	os << indent << "final_scaled_extension_gap_count = " << final_scaled_extension_gap_count << endl;
	os << indent << "final_scaled_extension_sequence_count = " << final_scaled_extension_sequence_count << endl;
	os << indent << "initial_scaled_gap_count = " << initial_scaled_gap_count << endl;
	os << indent << "initial_scaled_sequence_count = " << initial_scaled_sequence_count << endl;
	os << indent << "initial_scaled_extension_gap_count = " << initial_scaled_extension_gap_count << endl;
	os << indent << "initial_scaled_extension_sequence_count = " << initial_scaled_extension_sequence_count << endl;
	os << indent << "scaled_gap_deletion_penalty = " << scaled_gap_deletion_penalty << endl;
	os << indent << "scaled_gap_deletion_penalty_opens = " << scaled_gap_deletion_penalty_opens << endl;
	os << indent << "final_gap_closing_penalty = " << final_gap_closing_penalty << endl;
	os << indent << "initial_gap_opening_penalty = " << initial_gap_opening_penalty << endl;
	
	os << indent << "initial_sequence_count = " << initial_sequence_count << endl;
	os << indent << "initial_gap_count = " << initial_gap_count << endl;
	os << indent << "initial_scaled_gap_deletion_penalty = " << initial_scaled_gap_deletion_penalty << endl;
	// os << indent << "initial_scaled_gap_deletion_penalty_opens = " << initial_scaled_gap_deletion_penalty_opens << endl;

	os << indent << "gap_opening_penalty = " << gap_opening_penalty << endl;
	os << indent << "gap_closing_penalty = " << gap_closing_penalty << endl;
	os << indent << "gap_deletion_penalty = " << gap_deletion_penalty << endl;
	os << indent << "is_initial = " << is_initial << endl;
	os << indent << "is_final = " << is_final << endl;
	
	// os << indent << " gap parameter weight = " << p_gap_parameter.sequence_weight << endl;
	
	os << indent << "Residue counts:" << endl;
	map<char,float>::const_iterator count;
	for (count = residue_counts.begin(); count != residue_counts.end(); count++) {
		os << indent << "  " << count->first << ": " << count->second << endl;
	}
	os << indent << "Residue score counts:" << endl;
	for (char r = 'A'; r <= 'Z'; r++) {
		os << indent << "  " << r << ": " << get_residue_score_count(r) << endl;
	}
	os << indent << "Sequence residues:" << endl;
	for (unsigned int i = 0; i < sequence_residues.size(); i++) {
		os << indent << "  sequence " << i << ", residue " << sequence_residues[i] << endl;
	}
	
	return os;
}

void Conservidue::initialize_members() {
	gap_opening_penalty.opening() = 0;
	gap_closing_penalty.closing() = 0;
	gap_deletion_penalty.deletion() = 0;
	parent_alignment = NULL;
	weighted_sequence_count = 0;
	weighted_gap_count = 0;

	scaled_extension_gap_count = 0;
	scaled_gap_count = 0;
	scaled_gap_open_count = 0;
	scaled_gap_close_count = 0;
	scaled_extension_left_sequence_count = weighted_sequence_count;
	scaled_left_sequence_count = weighted_sequence_count;
	scaled_extension_right_sequence_count = weighted_sequence_count;
	scaled_right_sequence_count = weighted_sequence_count;
	
	scaled_extension_gap_open_count = 0;
	scaled_extension_gap_close_count = 0;
	final_scaled_gap_count = 0;
	final_scaled_sequence_count = 0;
	final_scaled_extension_gap_count = 0;
	final_scaled_extension_sequence_count = 0;
	initial_scaled_gap_count = 0;
	initial_scaled_sequence_count = 0;
	initial_scaled_extension_gap_count = 0;
	initial_scaled_extension_sequence_count = 0;
	scaled_gap_deletion_penalty.deletion() = 0;
	scaled_gap_deletion_penalty_opens.deletion() = 0;
	final_gap_closing_penalty.closing() = 0;
	initial_gap_opening_penalty.closing() = 0;

	initial_sequence_count = 0;
	initial_gap_count = 0;
	initial_scaled_gap_deletion_penalty.deletion() = 0;
	// initial_scaled_gap_deletion_penalty_opens.deletion() = 0;

	residue_score_counts.assign(26, 0);		
	array_sequence_index = 0;
	is_initial = false;
	is_final = false;
}

Conservidue::Conservidue() : 
parent_alignment(NULL),
weighted_sequence_count(0.0),
weighted_gap_count(0.0),

scaled_extension_gap_count(0),
scaled_gap_count(0),
scaled_gap_open_count(0),
scaled_gap_close_count(0),
scaled_extension_left_sequence_count(weighted_sequence_count),
scaled_left_sequence_count(weighted_sequence_count),
scaled_extension_right_sequence_count(weighted_sequence_count),
scaled_right_sequence_count(weighted_sequence_count),

scaled_extension_gap_open_count(0),
scaled_extension_gap_close_count(0),
final_scaled_gap_count(0),
final_scaled_sequence_count(0),
final_scaled_extension_gap_count(0),
final_scaled_extension_sequence_count(0),
initial_scaled_gap_count(0),
initial_scaled_sequence_count(0),
initial_scaled_extension_gap_count(0),
initial_scaled_extension_sequence_count(0),
initial_sequence_count(0),
initial_gap_count(0),
initial_scaled_gap_deletion_penalty(DELETION_SCORE, 0),
// initial_scaled_gap_deletion_penalty_opens(DELETION_SCORE, 0),
scaled_gap_deletion_penalty(DELETION_SCORE, 0),
scaled_gap_deletion_penalty_opens(DELETION_SCORE, 0),
final_gap_closing_penalty(CLOSING_SCORE, 0),
initial_gap_opening_penalty(OPENING_SCORE, 0),

gap_opening_penalty(OPENING_SCORE, 0),
gap_closing_penalty(CLOSING_SCORE, 0),
gap_deletion_penalty(DELETION_SCORE, 0),
p_gap_parameter(protein_gap_model, weighted_sequence_count),
residue_score_counts(26, 0),
array_sequence_index(0),
is_initial(false),
is_final(false)
{
	// initialize_members();
} // Simple constructor	

Conservidue::Conservidue(const Residue & residue) :
parent_alignment(NULL),
weighted_sequence_count(residue.sequence_pointer()->get_weight()),
weighted_gap_count(0.0),

scaled_extension_gap_count(0),
scaled_gap_count(0),
scaled_gap_open_count(0),
scaled_gap_close_count(0),
scaled_extension_left_sequence_count(weighted_sequence_count),
scaled_left_sequence_count(weighted_sequence_count),
scaled_extension_right_sequence_count(weighted_sequence_count),
scaled_right_sequence_count(weighted_sequence_count),

scaled_extension_gap_open_count(0),
scaled_extension_gap_close_count(0),
final_scaled_gap_count(0),
final_scaled_sequence_count(0),
final_scaled_extension_gap_count(0),
final_scaled_extension_sequence_count(0),
initial_scaled_gap_count(0),
initial_scaled_sequence_count(0),
initial_scaled_extension_gap_count(0),
initial_scaled_extension_sequence_count(0),

initial_sequence_count(0),
initial_gap_count(0),
initial_scaled_gap_deletion_penalty(DELETION_SCORE, 0),
// initial_scaled_gap_deletion_penalty_opens(DELETION_SCORE, 0),

scaled_gap_deletion_penalty(DELETION_SCORE, 0),
scaled_gap_deletion_penalty_opens(DELETION_SCORE, 0),
final_gap_closing_penalty(CLOSING_SCORE, 0),
initial_gap_opening_penalty(OPENING_SCORE, 0),

gap_opening_penalty(OPENING_SCORE, 0),
gap_closing_penalty(CLOSING_SCORE, 0),
gap_deletion_penalty(DELETION_SCORE, 0),
p_gap_parameter(protein_gap_model, weighted_sequence_count),
residue_score_counts(26, 0),
array_sequence_index(0),
is_initial(false),
is_final(false)
{
	const BioSequence * sequence = residue.sequence_pointer();
	
	residue_counts[residue.one_letter_code()] = weighted_sequence_count;

	gap_opening_penalty = residue.p_gap_opening_penalty * weighted_sequence_count;
	gap_closing_penalty = residue.p_gap_closing_penalty * weighted_sequence_count;
	gap_deletion_penalty = residue.p_gap_deletion_penalty * weighted_sequence_count;
	scaled_gap_deletion_penalty = residue.p_gap_deletion_penalty * weighted_sequence_count;
	
	// residue_score_count for efficient match score computation
	for (char res2 = 'A'; res2 <= 'Z'; res2 ++) {
		char res1 = residue.one_letter_code();
		float bit_score = blosum62.get_score(res1,res2);
		residue_score_count(res2) = bit_score * weighted_sequence_count;
	}
	
	int sequence_number = sequence->alignment_sequence_index;
	int residue_index = residue.sequence_residue_index;
	if (sequence_number == 0) { // create new sequence_residues mapping
		sequence_residues.clear();
		sequence_residues.push_back(residue_index);
	}
	else {
		sequence_residues.assign(sequence_number, -1); // default residue number -1 for other sequences
		sequence_residues[sequence_number] = residue_index;
	}
	
	// Efficiency members for scoring gaps in profile alignments
	int final_residue = sequence->length() - 1;
	if (residue_index == 0) { // first residue in sequence
		scaled_extension_left_sequence_count *= sequence->left_gap_extension_factor;
		scaled_left_sequence_count *= sequence->left_gap_factor;
		// This is the one case where initial_gap_opening_penalty is non-zero
		initial_gap_opening_penalty = gap_opening_penalty;
		initial_scaled_sequence_count = sequence->left_gap_factor * weighted_sequence_count;
		initial_sequence_count = weighted_sequence_count;
		initial_scaled_extension_sequence_count = sequence->left_gap_extension_factor * weighted_sequence_count;
		// initial_scaled_gap_deletion_penalty = residue.p_gap_deletion_penalty * weighted_sequence_count * sequence->left_gap_factor;
	}
	if (residue_index == final_residue) { // final residue in sequence
		scaled_extension_right_sequence_count *= sequence->right_gap_extension_factor;
		scaled_right_sequence_count *= sequence->right_gap_factor;
		scaled_gap_deletion_penalty *= sequence->right_gap_factor;
		// This is the one case where final_gap_closing_penalty is non-zero
		final_gap_closing_penalty = gap_closing_penalty;
		final_scaled_sequence_count = sequence->right_gap_factor * weighted_sequence_count;
		final_scaled_extension_sequence_count = sequence->right_gap_extension_factor * weighted_sequence_count;
	}

}

