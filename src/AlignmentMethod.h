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
