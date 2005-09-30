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
