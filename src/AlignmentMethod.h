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


#ifndef __ALIGNMENT_METHOD_H__
#define __ALIGNMENT_METHOD_H__

// Header for alignment methods used to align objects defined in SequenceAlignment.h

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
		align_global(true), use_divergence(NO_DIVERGENCE), branch_alignment(false)
	{}
};

extern AlignmentGranularity default_granularity;

// declaration needed only so it can be called friend to SequenceAlignment
class ConservidueAlignment; // single cell of dynamic programming table

#include "SequenceAlignment.h"

// Basic sequence alignment rountine, for a single step of progressive alignment
// This is a wrapper for SequenceAlignment::align()
SequenceAlignment align_profiles(
	const SequenceAlignment & seq1, 
	const SequenceAlignment & seq2,
	AlignmentGranularity granularity = default_granularity
								 );

#endif
