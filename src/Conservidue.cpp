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

// $Id$
// $Header$
// $Log$
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

Conservidue Conservidue::indel_conservidue() const {
	Conservidue answer;
	
	const Conservidue * before = this;
	const Conservidue * after = NULL;
	if (!before->is_final) {
		after = &(before->parent_alignment->operator[](before->array_sequence_index + 1));
	}
	
	if ((before->is_initial) && (after == NULL)) throw 99;
	
	const SequenceAlignment * sequence_alignment = NULL;
	
	if (before->is_initial) { // Left end gap
		sequence_alignment = after->parent_alignment;
		answer.weighted_left_end_gap_count = sequence_alignment->weighted_sequence_count;
		answer.weighted_right_end_gap_count = 0;
		answer.weighted_internal_gap_count = 0;
	}
	else if (after == NULL) { // Right end gap
		sequence_alignment = before->parent_alignment;
		answer.weighted_left_end_gap_count = 0;
		answer.weighted_right_end_gap_count = sequence_alignment->weighted_sequence_count;
		answer.weighted_internal_gap_count = 0;
	}
	else {
		sequence_alignment = before->parent_alignment;
		answer.weighted_left_end_gap_count = after->weighted_left_end_gap_count + after->weighted_left_end_sequence_count;
		answer.weighted_right_end_gap_count = before->weighted_right_end_gap_count + before->weighted_right_end_sequence_count;
		// ALERT assumes that all sequences are at least 2 residues long
		answer.weighted_internal_gap_count = before->weighted_internal_sequence_count + before->weighted_left_end_sequence_count;
	}
	
	// TODO - gap open count depends upon whether before conservidue is EXACTLY before
	answer.gap_open_count = 0;
	// Pure indels can not be gap close positions
	answer.gap_close_count = 0;
	
	answer.weighted_sequence_count = 0;
	answer.sequence_residues.assign(sequence_alignment->sequences.size(), -1);
	
	return answer;
}

// Simple sum of pairs formula
// A(x,y) = match score = SUM(over i,j) [res(i,x) * res(j,y) * Score(i,j)]
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

// C(x,y) = match gap extension score = EP(y) * gap_count(x) + EP(x) * gap_count(y)
AlignmentScore Conservidue::match_gap_extension_score(const Conservidue & c2, unsigned int gseg) const {
	const Conservidue & c1 = *this;	
	AlignmentScore answer;

	const SequenceAlignment & seq1 = *(c1.parent_alignment);
	const SequenceAlignment & seq2 = *(c2.parent_alignment);

	// TODO - precompute these factors
	answer += c1.gap_extension_penalty(gseg) *
		(c2.weighted_internal_gap_count +
		 c2.weighted_left_end_gap_count * seq2.left_gap_extension_factor +
		 c2.weighted_right_end_gap_count * seq2.right_gap_extension_factor) +
		c2.gap_extension_penalty(gseg) *
		(c1.weighted_internal_gap_count +
		 c1.weighted_left_end_gap_count * seq1.left_gap_extension_factor +
		 c1.weighted_right_end_gap_count * seq1.right_gap_extension_factor);		
	
	return answer;
}

Conservidue::Conservidue(const Residue & residue) :
parent_alignment(NULL),
weighted_sequence_count(residue.sequence_pointer()->get_weight()),
weighted_internal_gap_count(0),
weighted_left_end_gap_count(0),
weighted_right_end_gap_count(0),
weighted_left_end_sequence_count(0),
weighted_right_end_sequence_count(0),
weighted_internal_sequence_count(weighted_sequence_count),
gap_open_count(0),			
gap_close_count(0),			
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
	
	// initialize_members();
	// weighted_sequence_count = sequence->weight;
	// weighted_internal_gap_count = 0.0;
	// weighted_left_end_gap_count = 0.0;
	// weighted_right_end_gap_count = 0.0;
	
	residue_counts[residue.one_letter_code()] = weighted_sequence_count;
	gap_opening_penalty = residue.p_gap_opening_penalty * weighted_sequence_count;
	gap_closing_penalty = residue.p_gap_closing_penalty * weighted_sequence_count;
	gap_deletion_penalty = residue.p_gap_deletion_penalty * weighted_sequence_count;

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
	answer.weighted_internal_gap_count += conservidue2.weighted_internal_gap_count;
	answer.weighted_left_end_gap_count += conservidue2.weighted_left_end_gap_count;
	answer.weighted_right_end_gap_count += conservidue2.weighted_right_end_gap_count;
	answer.weighted_left_end_sequence_count += conservidue2.weighted_left_end_sequence_count;
	answer.weighted_right_end_sequence_count += conservidue2.weighted_right_end_sequence_count;
	answer.weighted_internal_sequence_count += conservidue2.weighted_internal_sequence_count;
	answer.gap_open_count += conservidue2.gap_open_count;
	answer.gap_close_count += conservidue2.gap_close_count;
	
	answer.gap_opening_penalty += conservidue2.gap_opening_penalty;
	answer.gap_closing_penalty += conservidue2.gap_closing_penalty;
	answer.p_gap_parameter = conservidue1.p_gap_parameter.combine(conservidue2.p_gap_parameter);
	answer.gap_deletion_penalty += conservidue2.gap_deletion_penalty;
	answer.parent_alignment = NULL;
	
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

void Conservidue::initialize_members() {
	gap_opening_penalty.opening() = 0;
	gap_closing_penalty.closing() = 0;
	gap_deletion_penalty.deletion() = 0;
	parent_alignment = NULL;
	weighted_sequence_count = 0;
	weighted_internal_gap_count = 0;
	weighted_left_end_gap_count = 0;
	weighted_right_end_gap_count = 0;
	weighted_left_end_sequence_count = 0;
	weighted_right_end_sequence_count = 0;
	weighted_internal_sequence_count = 0;
	gap_open_count = 0;
	gap_close_count = 0;
	residue_score_counts.assign(26, 0),
	array_sequence_index = 0;
	is_initial = false;
	is_final = false;
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
	os << indent << "weighted_internal_gap_count = " << weighted_internal_gap_count << endl;
	os << indent << "weighted_left_end_gap_count = " << weighted_left_end_gap_count << endl;
	os << indent << "weighted_right_end_gap_count = " << weighted_right_end_gap_count << endl;
	os << indent << "weighted_left_end_sequence_count = " << weighted_left_end_sequence_count << endl;
	os << indent << "weighted_right_end_sequence_count = " << weighted_right_end_sequence_count << endl;
	os << indent << "weighted_internal_sequence_count = " << weighted_internal_sequence_count << endl;
	os << indent << "gap_open_count = " << gap_open_count << endl;
	os << indent << "gap_close_count = " << gap_close_count << endl;
	
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

Conservidue::Conservidue() : 
parent_alignment(NULL),
weighted_sequence_count(0.0),
weighted_internal_gap_count(0),
weighted_left_end_gap_count(0),
weighted_right_end_gap_count(0),
weighted_left_end_sequence_count(0),
weighted_right_end_sequence_count(0),
weighted_internal_sequence_count(0),
gap_open_count(0),
gap_close_count(0),
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
