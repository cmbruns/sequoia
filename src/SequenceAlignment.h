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

using namespace std; // So STL stuff will work in the traditional way

#define DEFAULT_LEFT_GAP_FACTOR 1.0
#define DEFAULT_RIGHT_GAP_FACTOR 1.0
#define DEFAULT_LEFT_GAP_EXTENSION_FACTOR 1.0
#define DEFAULT_RIGHT_GAP_EXTENSION_FACTOR 1.0

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
	Conservidue * predecessor_conservidue;
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
	vector<const Residue *> residues;
	map<const BioSequence *, const Residue *> sequence_residues; // ALERT - assumes that each conservidue has at most one residue from each sequence

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
	// Given a sequence, return the residue that is in this conservidue (for printing)
	const Residue * sequence_residue(const BioSequence * seq_ptr) const;
	map<char, float> residue_counts; // weighted by sequence weights
	bool is_initial; // Conservidue has no predecessors, fake Conservidue to indicate start state
	bool is_final; // Conservidue has no successors, fake Conservidue to indicate stop state

	ostream & print(ostream & os, unsigned int indent_size) const {
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
	// simple string output routine
	friend ostream & operator<<(ostream & os, const SequenceAlignment & s);
protected:
	vector<Conservidue> conservidues; // ALERT - every Conservidue's precursors must appear before the Conservidue itself
	vector<Conservidue *> begins; 
	vector<Conservidue *> ends;
	Macromolecule macromolecule; // of PROTEIN, DNA, RNA, XNA/NUCLEOTIDE, MRNA/CDNA, GENOMIC
	vector<const BioSequence *> sequences;
	
	// Scale scores of end gaps by this amount, but not extensions
	double left_gap_factor; // coefficient for left end gaps, other than extensions
	double right_gap_factor; // coefficient for right end gaps, other than extensions
	// Additional scale factor applies only to extension portion of gap penalty
	double left_gap_extension_factor;
	double right_gap_extension_factor;
	
public:
	SequenceAlignment() : 
		left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
		right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
		left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
		right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR)
	{} // simple constructor

	SequenceAlignment(const BioSequence & seq);
	
	int length() const {return conservidues.size();}
	const Conservidue & operator[](int i) const {return conservidues.at(i);}
};

istream & operator>>(istream & is, SequenceAlignment & alignment) {
	// Read a series of sequences, and paste them together into an alignment
	while (is.good()) {
		BioSequence sequence;
		is >> sequence;
	}
	return is;
}

#endif
