/*
 *  SequenceAlignment.h
 *  sequoia4
 *
 *  Basic classes for Sequences
 *  1) Sequence alignment - one sequence or an alignment of sequences
 *  2) Conservidue - one residue or an alignment of residues
 *
 *  Sequence alignment METHODS are actually in ConservidueAlignment.h
 *
 *  Created by Christopher Bruns on Mon May 03 2004.
 *  Copyright (c) 2004  All rights reserved.
 *
 */

#ifndef __SEQUENCE_ALIGNMENT_H__
#define __SEQUENCE_ALIGNMENT_H__

#include <vector>
#include <map>
#include <iostream>
#include "BioSequence.h"
#include "AlignmentMethod.h"
#include "Exceptions.h"

using namespace std; // So STL stuff will work in the traditional way

#define DEFAULT_LEFT_GAP_FACTOR 0.5
#define DEFAULT_RIGHT_GAP_FACTOR 0.5
#define DEFAULT_LEFT_GAP_EXTENSION_FACTOR 0.0
#define DEFAULT_RIGHT_GAP_EXTENSION_FACTOR 0.0

// Alignment space determines the alignment method.
// For protein sequences, it should always be protein
//  but for nucleotide sequences, it could be protein or nucleic.
// Hybrid alignments could be created using two versions of a DNA sequence,
//  one structured for protein, one for nucleic, with appropriate transitions
//  between them
enum AlignmentSpace {
	PROTEIN_ALIGNMENT,
	NUCLEIC_ALIGNMENT
};

class Conservidue; // short declaration

// Class to capture the transition between one residue and its predecessors
class ConserviduePredecessor {
	friend class Conservidue;
	friend class SequenceAlignment;
	friend class ConservidueAlignment;
protected:
	unsigned int predecessor_conservidue;
	double transition_score; // log2(probability of taking this path), 0.0 for a normal sequence
};

// A single residue or a set of aligned residues in a sequence alignment.
// Corresponds to one column of a multiple sequence alignment
// TODO - make this conservidue contain multiple residues
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
	//		   O       C
	//          EEEEEEEE

	double gap_opening_penalty; // -log2 probability of loop beginning with this Conservidue
	double gap_closing_penalty; // -log2 probability of loop before this Conservidue
	double gap_extension_penalty; // -log2 probability of this residue adding to a loop
	double gap_deletion_penalty; // -log2 probability of gap after this Conservidue
	
	SequenceAlignment * parent_alignment;
	
	vector<ConserviduePredecessor> predecessors;

	float weighted_sequence_count; // sum over all sequence weights
	// vector<const Residue *> residues; // do we need this, with sequence residues?

	// There is a public accessor for this
	// array index is sequence number in alignment
	// array value is residue number in sequence
	// value -1 means no such sequence in conservidue
	vector<int> sequence_residues; // ALERT - assumes that each conservidue has at most one residue from each sequence

	void initialize_members() {
		gap_opening_penalty = -1.1;
		gap_closing_penalty = -1.1;
		gap_extension_penalty = -0.23;
		gap_deletion_penalty = -1.0;
		parent_alignment = NULL;
		weighted_sequence_count = 0;
		array_sequence_index = 0;
		is_initial = false;
		is_final = false;
	}
public:
	int array_sequence_index; // actual order of Conservidue in parent sequence
	map<char, float> residue_counts; // weighted by sequence weights
	bool is_initial; // Conservidue has no predecessors, fake Conservidue to indicate start state
	bool is_final; // Conservidue has no successors, fake Conservidue to indicate stop state

	// Given a sequence, return the residue that is in this conservidue (for printing)
	const Residue * sequence_residue(unsigned int sequence_number) const;

	Conservidue combine_conservidues(const Conservidue & conservidue2) const; // for aligned conservidues
	
	ostream & print_debug(ostream & os = cout, unsigned int indent_size = 0) const {
		string indent = "";
		for(unsigned int i=0;i<indent_size;i++)indent += " ";

		os << indent << "conservidue pointer = " << this << endl;
		os << indent << "parent_alignment = " << parent_alignment << endl;
		os << indent << "array_sequence_index = " << array_sequence_index << endl;
		os << indent << "weighted_sequence_count = " << weighted_sequence_count << endl;
		os << indent << "gap_opening_penalty = " << gap_opening_penalty << endl;
		os << indent << "gap_closing_penalty = " << gap_closing_penalty << endl;
		os << indent << "gap_extension_penalty = " << gap_extension_penalty << endl;
		os << indent << "gap_deletion_penalty = " << gap_deletion_penalty << endl;
		os << indent << "is_initial = " << is_initial << endl;
		os << indent << "is_final = " << is_final << endl;

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
	
	Conservidue(const Residue & residue); 
	Conservidue() {
		initialize_members();
	} // Simple constructor	
};

// SequenceAlignment is both the input and output of progressive multiple sequence
// alignmnent.  It is a directed graph of conservidues.  An individual sequence
// is cast as a SequenceAlignment before it get aligned to other stuff.
class SequenceAlignment {
	friend class ConservidueAlignment;
	friend class Conservidue;
	// simple string output routine
	friend ostream & operator<<(ostream & os, const SequenceAlignment & s);
protected:
	// Scale scores of end gaps by this amount, but not extensions
	double left_gap_factor; // coefficient for left end gaps, other than extensions
	double right_gap_factor; // coefficient for right end gaps, other than extensions
	// Additional scale factor applies only to extension portion of gap penalty
	double left_gap_extension_factor;
	double right_gap_extension_factor;	

	vector<Conservidue> conservidues; // ALERT - every Conservidue's precursors must appear before the Conservidue itself
	vector<unsigned int> begins; 
	vector<unsigned int> ends;
	vector<BioSequence> sequences;
	
	double pair_alignment_score; // Score resulting from most recent alignment
								 // as opposed to sum of pairs score

public:
	int length() const {return conservidues.size();}
	const Conservidue & operator[](int i) const {return conservidues.at(i);}
	
	// Core profile alignment routine
	// See AlignmentMethod.cpp for implementation
	SequenceAlignment SequenceAlignment::align(const SequenceAlignment & seq2) const;

	void add_sequence(const BioSequence & sequence) { // does not automatically update conservidues!!
		sequences.push_back(sequence);
		int sequence_index = sequences.size() - 1;
		sequences.back().alignment_sequence_index = sequence_index;
		sequences.back().parent_alignment = this;
	}
	
	void add_conservidue(const Conservidue & c, unsigned int sequence_index_offset = 0) {
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
		pred.transition_score = 0;
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

	ostream & print_debug(ostream & os = cout, unsigned int indent_size = 0) const {
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
	
	SequenceAlignment() : 
		left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
		right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
		left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
		right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
		pair_alignment_score(0)
	{} // simple constructor
	SequenceAlignment(const BioSequence & seq);

	// need assignment operator and copy constructor so that conservidue.parent_alignment pointers are properly set
	SequenceAlignment(const SequenceAlignment & alignment2) {
		*this = alignment2;
	}
	
	SequenceAlignment & operator=(const SequenceAlignment & alignment2) {
		if (this == &alignment2) return *this;
		begins = alignment2.begins;
		ends = alignment2.ends;
		left_gap_factor = alignment2.left_gap_factor;
		right_gap_factor = alignment2.right_gap_factor;
		left_gap_extension_factor = alignment2.left_gap_extension_factor;
		right_gap_extension_factor = alignment2.right_gap_extension_factor;
		pair_alignment_score = alignment2.pair_alignment_score;

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
};

istream & operator>>(istream & is, SequenceAlignment & alignment);

#endif
