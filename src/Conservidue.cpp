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

// $Id$
// $Header$
// $Log$
// Revision 1.1  2004/05/29 22:40:21  cmbruns
// Created separate files for Conservidue, separate from SequenceAlignment
//

#include "Conservidue.h"
#include "SequenceAlignment.h"

// A single residue or a set of aligned residues in a sequence alignment.
// Corresponds to one column of a multiple sequence alignment

// ********* Conservidue methods *************

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
array_sequence_index(0),
is_initial(false),
is_final(false)
{
	// initialize_members();
} // Simple constructor	
