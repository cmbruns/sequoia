#ifndef __ALIGNMENT_METHOD_H__
#define __ALIGNMENT_METHOD_H__

// Header for alignment methods used to align objects defined in SequenceAlignment.h

enum AlignmentGranularity {
	ALIGN_GLOBAL, // Needleman-Wunsch style
	ALIGN_LOCAL // Smith-Waterman style
};

// declaration needed only so it can be called friend to SequenceAlignment
class ConservidueAlignment; // single cell of dynamic programming table

#include "SequenceAlignment.h"

// Basic sequence alignment rountine, for a single step of progressive alignment
// Should this be a member function of some class? -> it is now, SequenceAlignment
// but the actual method is defined in AlignmentMethod.cpp
SequenceAlignment align_profiles(
	const SequenceAlignment & seq1, 
	const SequenceAlignment & seq2,
	AlignmentGranularity granularity = ALIGN_GLOBAL
								 );

#endif
