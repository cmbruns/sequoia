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

// ConservidueAlignment object corresponds to one cell of a dynamic programing matrix

#ifndef __CONSERVIDUE_ALIGNMENT_H__
#define __CONSERVIDUE_ALIGNMENT_H__

#include <iostream>
#include <vector>
#include "AlignmentScore.h"

class DPMatrix;

// How to model regions of no relationship between two sequences
enum DivergenceType {
	NO_DIVERGENCE, // No special model
	PAIR_DIVERGENCE, // Explicit alignment step for pairs of unalignable residues
	GAP_DIVERGENCE // Model divergence using free pass between inserts and deletes
};

class AlignmentGranularity {
public:
	bool align_global; // true if Needleman-Wunsch, false if Smith-Waterman
	DivergenceType use_divergence; // Permit unaligned pairs?  LOBAL
	bool branch_alignment; // like POA
	
	// Set default alignment method here
	AlignmentGranularity() :
		align_global(true), use_divergence(GAP_DIVERGENCE), branch_alignment(false)
	{}
};

extern AlignmentGranularity default_granularity;

extern double BAD_SCORE;

// AlignmentStep is contained in the ConservidueAlignment class
// Each AlignmentStep represents one step in a traceback path
// of the dynamic programming table
// There is a separate AlignmentStep for each potential match, deletion, insertion,
// (and others?) state in the alignment path.
class AlignmentStep {
	friend class ConservidueAlignment;
protected:
public:
	const AlignmentStep * traceback_pointer;
	AlignmentScore path_score;
	const ConservidueAlignment * parent_cell;
	bool is_match; // true if this step actually aligns two residues, i.e. match state
	
	ostream & print_debug(ostream & os, unsigned int indent_size) const;

	AlignmentStep();
	AlignmentStep(const ConservidueAlignment * parent);
};

// ConservidueAlignment is a single cell in a dynamic programming
// sequence alignment table.  It relates one Residue/profile in the
// first sequence, to one in the second.  For example, one ConservidueAlignment
// might represent the prospect of alignment residue VAL 123 of sperm whale
// myoglobin to residue LYS 17 if human hemoglobin.
class ConservidueAlignment {
protected:
public:
	vector<AlignmentStep> deletion; // delete is a C++ keyword..., gap in sequence 1, loop in 2, "E" by Gusfield
	vector<AlignmentStep> insertion; // insertion path step(s), gap in sequence 2, loop in 1, "F" by Gusfield
	vector<AlignmentStep> divergence; // consume both Conservidues, but do not align them (for POA)
	AlignmentStep match; // no gap, alignment of sequence1 to sequence 2
	AlignmentStep * best; // highest scoring path step
	
	DPMatrix * parent_table;
	const Conservidue * conservidue1; // Conservidue/profile from first sequence
	const Conservidue * conservidue2; // Conservidue/profile from second sequence
	AlignmentScore conservidue_score; // Log odds score for aligning these two Conservidues
	
	void alignment_recurrence(AlignmentGranularity granularity); // Dynamic programming recurrence
	
	void initialize_upper_left(AlignmentGranularity granularity);
	void initialize_left_column(AlignmentGranularity granularity);
	void initialize_top_row(AlignmentGranularity granularity);
	void compute_pair_score();
	void assign_best_match_state(AlignmentGranularity granularity);
	void assign_best_deletion_state(AlignmentGranularity granularity);
	void assign_best_insertion_state(AlignmentGranularity granularity);
	void assign_best_divergence_state(AlignmentGranularity granularity);
	void assign_best_best_state(AlignmentGranularity granularity);
	
	const ConservidueAlignment * get_cell(const Conservidue * c1, const Conservidue * c2) const;	
	ostream & print_debug(ostream & os, unsigned int indent_size) const;
	
	ConservidueAlignment & operator=(const ConservidueAlignment & c2);
	ConservidueAlignment();	
	ConservidueAlignment(const ConservidueAlignment & c2);
};

#endif
