#ifndef __ALIGNMENT_METHOD_H__
#define __ALIGNMENT_METHOD_H__

// Header for alignment methods used to align objects defined in SequenceAlignment.h

#include "SequenceAlignment.h"

// Basic sequence alignment rountine, for a single step of progressive alignment
// Should this be a member function of some class?
SequenceAlignment align(const SequenceAlignment & seq1, const SequenceAlignment & seq2);

class ConservidueAlignment; // single cell of dynamic programming table

#endif
