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
	// simple string output routine
	friend ostream & operator<<(ostream & os, const Conservidue & c);
protected:
	double gap_open_penalty; // -log2 probability of gap after this Conservidue
	double gap_close_penalty; // -log2 probability of gap before this Conservidue
	int begin_distance; // distance to nearest begin Conservidue (for end gaps free only)
	bool is_begin; // Conservidue has no predecessors, fake Conservidue to indicate start state
	bool is_final; // Conservidue has no successors, fake Conservidue to indicate stop state
	int array_sequence_index; // actual order of Conservidue in parent sequence
	vector<ConserviduePredecessor> predecessors;

	float weighted_sequence_count; // sum over all sequence weights
	vector<const Residue *> residues;
	map<const BioSequence *, const Residue *> sequence_residues; // ALERT - assumes that each conservidue has at most one residue from each sequence

public:
	Conservidue(const Residue & residue); 
	Conservidue() {} // Simple constructor
	
	// Given a sequence, return the residue that is in this conservidue (for printing)
	const Residue * sequence_residue(const BioSequence * seq_ptr) const;
	map<char, float> residue_counts; // weighted by sequence weights
};

// SequenceAlignment is both the input and output of progressive multiple sequence
// alignmnent.  It is a directed graph of conservidues.  An individual sequence
// is cast as a SequenceAlignment before it get aligned to other stuff.
class SequenceAlignment {
	// simple string output routine
	friend ostream & operator<<(ostream & os, const SequenceAlignment & s);
protected:
	vector<Conservidue> conservidues; // ALERT - every Conservidue's precursors must appear before the Conservidue itself
	vector<Conservidue *> begins; 
	vector<Conservidue *> ends;
	Macromolecule macromolecule; // of PROTEIN, DNA, RNA, XNA/NUCLEOTIDE, MRNA/CDNA, GENOMIC
	vector<const BioSequence *> sequences;
public:
	SequenceAlignment() {} // simple constructor

	SequenceAlignment(const BioSequence & seq);
	
	int length() const {return conservidues.size();}
	const Conservidue & operator[](int i) const {return conservidues.at(i);}
};

#endif
