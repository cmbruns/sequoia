/*
 *  SequenceAlignment.h
 *  sequoia4
 *
 *  Basic classes for Sequences
 *  1) Sequence alignment - one sequence or an alignment of sequences
 *  2) Conservidue - one residue or an alignment of residues
 *
 *  Sequence alignment METHODS are actually in AlignmentMethod.h
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
#include <fstream>
#include "BioSequence.h"
#include "AlignmentMethod.h"
#include "Exceptions.h"
#include "GapModel.h"

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
	friend class ConservidueAlignment;
protected:
	unsigned int predecessor_conservidue;
	AlignmentScore transition_score; // log2(probability of taking this path), 0.0 for a normal sequence
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
	
	void initialize_members();

public:
	// Public variables
	int array_sequence_index; // actual order of Conservidue in parent sequence
	bool is_initial; // Conservidue has no predecessors, fake Conservidue to indicate start state
	bool is_final; // Conservidue has no successors, fake Conservidue to indicate stop state
	map<char, float> residue_counts; // weighted by sequence weights

	// Member functions
	Conservidue combine_conservidues(const Conservidue & conservidue2) const; // for aligned conservidues
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

// SequenceAlignment is both the input and output of progressive multiple sequence
// alignmnent.  It is a directed graph of conservidues.  An individual sequence
// is cast as a SequenceAlignment before it get aligned to other stuff.
class SequenceAlignment {
	friend class ConservidueAlignment;
	friend class Conservidue;
protected:
	// Scale scores of end gaps by this amount, but not extensions
	double left_gap_factor; // coefficient for left end gaps, other than extensions
	double right_gap_factor; // coefficient for right end gaps, other than extensions
	// Additional scale factor applies only to extension portion of gap penalty
	double left_gap_extension_factor;
	double right_gap_extension_factor;	

	double weighted_sequence_count;
	
	vector<Conservidue> conservidues; // ALERT - every Conservidue's precursors must appear before the Conservidue itself
	vector<unsigned int> begins; 
	vector<unsigned int> ends;
	vector<BioSequence> sequences;
	
	AlignmentScore pair_alignment_score; // Score resulting from most recent alignment
								 // as opposed to sum of pairs score

public:
	void add_conservidue(const Conservidue & c, unsigned int sequence_index_offset = 0);
	void add_sequence(const BioSequence & sequence);
	void add_sequence_automatic(const BioSequence & sequence);	
	// Core profile alignment routine
	// See AlignmentMethod.cpp for implementation
	SequenceAlignment align(const SequenceAlignment & seq2,
							AlignmentGranularity granularity = default_granularity) const;
	AlignmentScore alignment_score() const {return pair_alignment_score;}
	unsigned int length() const; 
	istream & SequenceAlignment::load_fasta(istream & is);
	void load_fasta_file(const char * file_name) {
		ifstream infile(file_name);
		if (infile == 0) {
			cerr << "*** ERROR *** : Unable to open file " << file_name << endl;
			throw NO_SUCH_FILE_EXCEPTION;
		}
		load_fasta(infile);
		infile.close();
	}
	const Conservidue & operator[](int i) const;
	ostream & print_debug(ostream & os = cout, unsigned int indent_size = 0) const;	
	ostream & SequenceAlignment::print_pretty(ostream & os = cout) const;
	float report_accuracy(const SequenceAlignment & true_alignment);
	vector<BioSequence> & sequence() {return sequences;}
	AlignmentScore sequence_pair_score(unsigned int seq_index1, unsigned int seq_index2) const; // TODO
	void set_gap_penalty(double penalty);
	void set_extension_penalty(const GapModel & gap_model);
	void set_weight(double w) {
		double factor = w / weighted_sequence_count; // Change everything by this factor
		weighted_sequence_count = w;
		for (unsigned int s = 0; s < sequences.size(); s++) {
			sequences[s].set_weight(sequences[s].get_weight() * factor);
		}
		for (unsigned int ci = 0; ci < conservidues.size(); ci ++) {
			Conservidue & c = conservidues[ci];

			c.weighted_sequence_count *= factor; // sum over all sequence weights with non-gaps

			c.weighted_internal_gap_count *= factor;
			c.weighted_left_end_gap_count *= factor;
			c.weighted_right_end_gap_count *= factor;
			c.weighted_left_end_sequence_count *= factor;  // count of initial residues
			c.weighted_right_end_sequence_count *= factor; // count of final residues
			c.weighted_internal_sequence_count *= factor; // count of internal residues			
			c.gap_opening_penalty *= factor; // -log2 probability of loop beginning with this Conservidue
			c.gap_closing_penalty *= factor; // -log2 probability of loop before this Conservidue
			c.gap_deletion_penalty *= factor; // -log2 probability of gap after this Conservidue

			c.p_gap_parameter.set_weight(c.weighted_sequence_count);
			
			map<char, float>::iterator res_count;
			for (res_count = c.residue_counts.begin(); 
				 res_count != c.residue_counts.end();
				 res_count ++) {
				res_count->second *= factor;
			}
		}
	}
	AlignmentScore sum_of_pairs_score() const;
	
	// need assignment operator and copy constructor so that conservidue.parent_alignment pointers are properly set
	SequenceAlignment & operator=(const SequenceAlignment & alignment2);
	SequenceAlignment();
	SequenceAlignment(const SequenceAlignment & alignment2);	
	SequenceAlignment(const BioSequence & seq);
};

istream & operator>>(istream & is, SequenceAlignment & alignment);
ostream & operator<<(ostream & os, const SequenceAlignment & s);

#endif
