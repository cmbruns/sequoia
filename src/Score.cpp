#include "Score.h"
#include "MutationMatrix.h"

float conservidue_score(const Conservidue & c1, const Conservidue & c2) {
	// Simple sum of pairs formula
	float answer = 0.0;
	
	map<char, float>::const_iterator res1;
	map<char, float>::const_iterator res2;
	for (res1 = c1.residue_counts.begin(); 
		 res1 != c1.residue_counts.end();
		 res1 ++) {
		for (res2 = c2.residue_counts.begin(); 
			 res2 != c2.residue_counts.end();
			 res2 ++) {
		}
		float sequence_count = res1->second * res2->second;
		float bit_score = blosum62.get_score(res1->first, res2->first);

		answer += sequence_count * bit_score;
	}
	
	return answer;
}

