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
// 

// $Id$
// $Header$
// $Log$
// Revision 1.7  2004/06/14 16:57:49  cmbruns
// Created update_begin_conservidue subroutine, but it should probably only be called in one place
// Improved logic for GAP_DIVERGENCE case in sequence_pair_score
// Created delta_sum_of_pairs_score for easier debugging in later progressive alignment steps
// Incorporated some new Conservidue members in set_weight() subroutine
//
// Revision 1.6  2004/06/04 19:31:10  cmbruns
// Updated GPL header
//
// Moved constructor methods to end of file
//
// Added initialize_from_pdb
//
// Moved all Conservidue functions to their own file Conservidue.h
//
// Enhanced pair score to work for local alignments, and to partially work for GAP_DIVERGENCE alignments.  TODO
//

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
// #include "AlignmentMethod.h"
#include "Exceptions.h"
// #include "GapModel.h"
#include "Conservidue.h"
#include "PDBEntry.h"
#include "ConservidueAlignment.h"

using namespace std; // So STL stuff will work in the traditional way

// TODO - eliminate sequences from SequenceAlignment
// TODO - eliminate end gap parameters from SequenceAlignment (leave in BioSequence)

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

// SequenceAlignment is both the input and output of progressive multiple sequence
// alignmnent.  It is a directed graph of conservidues.  An individual sequence
// is cast as a SequenceAlignment before it get aligned to other stuff.
class SequenceAlignment {
	friend class ConservidueAlignment;
	friend class Conservidue;
protected:
	// TODO - eliminate these gap factors, and let them propagate from BioSequence
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

	// I keep changing my mind about the contents of the begin conservidue, so make a function to take care of it
	void update_begin_conservidue();
	
	SequenceAlignment & initialize_from_pdb_protein(const PDBChain & chain);
	// void add_pdb_protein(PDBChain chain);
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
			throw NO_SUCH_FILE_EXCEPTION();
		}
		load_fasta(infile);
		infile.close();
	}
	// istream & SequenceAlignment::load_pdb(istream & is);
	void load_pdb_file(const char * file_name) {
		PDBEntry pdb_entry;
		pdb_entry.load_pdb_file(file_name);
		initialize_from_pdb_protein(pdb_entry.get_first_chain());
	}
	const Conservidue & operator[](int i) const;
	ostream & print_debug(ostream & os = cout, unsigned int indent_size = 0) const;	
	ostream & SequenceAlignment::print_pretty(ostream & os = cout) const;
	float report_accuracy(const SequenceAlignment & true_alignment);
	vector<BioSequence> & sequence() {return sequences;}
	AlignmentScore sequence_pair_score(unsigned int seq_index1, unsigned int seq_index2, AlignmentGranularity granularity = default_granularity) const; // TODO
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

			c.scaled_extension_gap_count *= factor;
			c.scaled_gap_count *= factor;
			c.scaled_gap_open_count *= factor;
			c.scaled_gap_close_count *= factor;
			c.scaled_extension_left_sequence_count *= factor;
			c.scaled_left_sequence_count *= factor;
			c.scaled_extension_right_sequence_count *= factor;
			c.scaled_right_sequence_count *= factor;
			c.scaled_extension_gap_open_count *= factor;
			c.scaled_extension_gap_close_count *= factor;
			c.scaled_gap_deletion_penalty *= factor;
			c.scaled_gap_deletion_penalty_opens *= factor;
			
			c.gap_opening_penalty *= factor; // -log2 probability of loop beginning with this Conservidue
			c.gap_closing_penalty *= factor; // -log2 probability of loop before this Conservidue
			c.gap_deletion_penalty *= factor; // -log2 probability of gap after this Conservidue

			c.p_gap_parameter.set_weight(c.weighted_sequence_count);
			
			// residue counts
			map<char, float>::iterator res_count;
			for (res_count = c.residue_counts.begin(); 
				 res_count != c.residue_counts.end();
				 res_count ++) {
				res_count->second *= factor;
			}
			
			// TODO - residue score counts
			vector<float>::iterator rs_count;
			for (rs_count = c.residue_score_counts.begin();
				 rs_count != c.residue_score_counts.end();
				 rs_count ++) {
				*rs_count *= factor;
			}
		}
	}
	AlignmentScore sum_of_pairs_score(AlignmentGranularity granularity = default_granularity) const;
	AlignmentScore delta_sum_of_pairs_score(unsigned int border_sequence, AlignmentGranularity granularity = default_granularity) const;
	SequenceAlignment & initialize_from_biosequence(const BioSequence & seq0);
	
	// need assignment operator and copy constructor so that conservidue.parent_alignment pointers are properly set
	SequenceAlignment & operator=(const SequenceAlignment & alignment2);
	SequenceAlignment();
	SequenceAlignment(const SequenceAlignment & alignment2);	
	SequenceAlignment(const BioSequence & seq);
	SequenceAlignment(const char * seq_string);
};

// Basic sequence alignment rountine, for a single step of progressive alignment
// This is a wrapper for SequenceAlignment::align()
SequenceAlignment align_profiles(
								 const SequenceAlignment & seq1, 
								 const SequenceAlignment & seq2,
								 AlignmentGranularity granularity = default_granularity
								 );

istream & operator>>(istream & is, SequenceAlignment & alignment);
ostream & operator<<(ostream & os, const SequenceAlignment & s);

#endif
