#ifndef __SCORE_H__
#define __SCORE_H__

#include "SequenceAlignment.h"

// return bit score of aligning two conservidues
float conservidue_pair_score(const Conservidue & c1, const Conservidue & c2);

#endif
