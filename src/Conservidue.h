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
// Revision 1.2  2004/06/04 19:06:10  cmbruns
// Updated GPL header
//
// Initial stubs for subroutines related to sum of pairs scoring in alignment
//
// Initial version of indel_conservidue subroutine, intended to help clean up the trace back code in AlignmentMethod.cpp
//
// Implemented residue_score_count member, to precompute scores for faster alignment
//
// Revision 1.1  2004/05/29 22:40:30  cmbruns
// Created separate files for Conservidue, separate from SequenceAlignment
//

#ifndef __CONSERVIDUE_H__
#define __CONSERVIDUE_H__

#include <map>
#include <vector>
#include "Score.h"
#include "GapModel.h"
#include "BioSequence.h"

// TODO - eliminate sequence_residues from Conservidue

// Class to capture the transition between one profile column and its predecessors
class ConserviduePredecessor {
	friend class Conservidue;
	friend class SequenceAlignment;
	friend class ConservidueAlignment;
protected:
		unsigned int predecessor_conservidue;
	AlignmentScore transition_score; // log2(probability of taking this path), 0.0 for a normal sequence
};

// A single residue or a set of aligned residues in a sequence alignment.
// Corresponds to one column of a multiple sequence alignment
class Conservidue {
	friend class SequenceAlignment;
	friend class ConservidueAlignment;
protected:
		
		//  Diagram of a typical alignment gap
		//  positions O, D, C, and E are positions where open, deletion, close,
		//  and extension penalties are computed
		//         D
		//			<-gap-->
		//  XXXXXXXX--------XXXXXXXXX gap sequence
		//  XXXXXXXXXXXXXXXXXXXXXXXXX loop sequence
		//          <-loop->
		//		    O      C
		//          EEEEEEEE
		
		SequenceAlignment * parent_alignment;
	
	vector<ConserviduePredecessor> predecessors;
	
	double weighted_sequence_count; // sum over all sequence weights with non-gaps
	double weighted_internal_gap_count;
	double weighted_left_end_gap_count;
	double weighted_right_end_gap_count;
	
	// left_end + right_end + internal = weighted_sequence_count
	double weighted_left_end_sequence_count;  // count of initial residues
	double weighted_right_end_sequence_count; // count of final residues
	double weighted_internal_sequence_count; // count of internal residues
	double gap_open_count; // sequences with gap here, amino acid before
	double gap_close_count; // sequences with gap here, amino acid after
	
	AlignmentScore gap_opening_penalty; // -log2 probability of loop beginning with this Conservidue
	AlignmentScore gap_closing_penalty; // -log2 probability of loop before this Conservidue
	AlignmentScore gap_deletion_penalty; // -log2 probability of gap after this Conservidue
	ResidueGapParameter p_gap_parameter;
	
	// vector<const Residue *> residues; // do we need this, with sequence residues?
	
	// There is a public accessor for this
	// array index is sequence number in alignment
	// array value is residue number in sequence
	// value -1 means no such sequence in conservidue
	vector<int> sequence_residues; // ALERT - assumes that each conservidue has at most one residue from each sequence
	vector<int> sequence_changes; // which sequences change their presence status vs.the previous residue
	
	map<char, float> residue_counts; // weighted by sequence weights
	vector<float> residue_score_counts; // precomputed score for each amino acid in alignment partner
	
	void initialize_members();
	
public:
	// Public variables
	int array_sequence_index; // actual order of Conservidue in parent sequence
	bool is_initial; // Conservidue has no predecessors, fake Conservidue to indicate start state
	bool is_final; // Conservidue has no successors, fake Conservidue to indicate stop state

	// Member functions

	// Populate an "empty" indel conservidue to be placed after the subject conservidue
	Conservidue indel_conservidue() const;
	
	// Precomputed score vector for efficient match score computation
	float get_residue_score_count(char amino_acid) const {
		int index = toupper(amino_acid - 'A');
		if (index < 0) return 0;
		if (index > ('Z' - 'A')) return 0;
		return residue_score_counts[index];
	}
	float & residue_score_count(char amino_acid) {
		int index = toupper(amino_acid - 'A');
		if (index < 0) throw 99;
		if (index > ('Z' - 'A')) throw 99;
		return residue_score_counts[index];
	}
	
	Conservidue combine_conservidues(const Conservidue & conservidue2) const; // for aligned conservidues

	// Support functions for alignment scores
	AlignmentScore match_score(const Conservidue & conservidue2) const;
	AlignmentScore match_gap_extension_score(const Conservidue & conservidue2, unsigned int gseg) const;

	// Piecewise linear gap model	
	AlignmentScore gap_extension_penalty(unsigned int gseg) const;
	AlignmentScore gap_open_offset(unsigned int gseg) const;
	AlignmentScore gap_close_offset(unsigned int gseg) const;
	unsigned int gap_segment_count() const {return p_gap_parameter.segment_count();}
	ostream & print_debug(ostream & os = cout, unsigned int indent_size = 0) const;
	// Given a sequence, return the residue that is in this conservidue (for printing)
	const Residue * sequence_residue(unsigned int sequence_number) const;
	
	// Constructors
	Conservidue(const Residue & residue); 
	Conservidue();
};


#endif

