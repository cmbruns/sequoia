#ifndef __ALIGNMENT_METHOD_H__
#define __ALIGNMENT_METHOD_H__

// Header for alignment methods used to align objects defined in SequenceAlignment.h

class SequenceAlignment;
class ConservidueAlignment; // single cell of dynamic programming table

#include "SequenceAlignment.h"

// Basic sequence alignment rountine, for a single step of progressive alignment
// Should this be a member function of some class?
SequenceAlignment align_profiles(const SequenceAlignment & seq1, const SequenceAlignment & seq2);

#endif
