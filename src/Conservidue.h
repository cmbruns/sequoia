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
#include "AlignmentScore.h"
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
	double weighted_gap_count; // TODO sum over all sequence weights with gaps
	
	// efficiency variables for sum of pairs
	// the variables with "scaled_" in the name are pre-adjusted for special end-gap weighting
	double scaled_extension_gap_count; // internal + left * factor + right * factor
	double scaled_gap_count; // needed?
	double scaled_gap_open_count; // needed?
	double scaled_gap_close_count;
	double scaled_extension_left_sequence_count; // also scaled by gap factor, like gap counts
	double scaled_left_sequence_count; // also scaled by gap factor, like gap counts
	double scaled_extension_right_sequence_count; // also scaled by gap factor, like gap counts
	double scaled_right_sequence_count; // also scaled by gap factor, like gap counts
	double scaled_extension_gap_open_count;
	double scaled_extension_gap_close_count;

	double final_scaled_gap_count;
	double final_scaled_sequence_count;
	double final_scaled_extension_gap_count;
	double final_scaled_extension_sequence_count;
	double initial_scaled_gap_count;
	double initial_scaled_sequence_count;
	double initial_scaled_extension_gap_count;
	double initial_scaled_extension_sequence_count;
	
	// For initial gap deletion
	double initial_sequence_count;
	double initial_gap_count; // TODO
	AlignmentScore initial_scaled_gap_deletion_penalty;
	// AlignmentScore initial_scaled_gap_deletion_penalty_opens; // TODO
	
	AlignmentScore scaled_gap_deletion_penalty;
	AlignmentScore scaled_gap_deletion_penalty_opens; // gap open positions only

	AlignmentScore final_gap_closing_penalty; // Only final residues of sequences
	AlignmentScore initial_gap_opening_penalty; // Only initial residues of sequences
	
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
	// vector<int> sequence_changes; // which sequences change their presence status vs.the previous residue
	
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
	Conservidue indel_conservidue(bool is_first_indel) const;
	
	// Precomputed score vector for efficient match score computation
	float get_residue_score_count(char amino_acid) const;
	float & residue_score_count(char amino_acid);
	
	Conservidue combine_conservidues(const Conservidue & conservidue2) const; // for aligned conservidues

	// Support functions for alignment scores
	AlignmentScore match_score(const Conservidue & conservidue2) const;
	AlignmentScore initial_gap_extension_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore match_gap_extension_score(const Conservidue & conservidue2, unsigned int gseg) const;
	AlignmentScore new_gap_extension_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore initial_gap_opening_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore match_gap_opening_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore new_gap_opening_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore new_gap_deletion_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore match_gap_closing_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore new_gap_closing_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore final_match_gap_closing_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore final_new_gap_closing_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore initial_match_gap_opening_score(const Conservidue & c2, unsigned int gseg) const;
	AlignmentScore initial_new_gap_opening_score(const Conservidue & c2, unsigned int gseg) const;

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

