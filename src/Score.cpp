#include "Score.h"
#include "MutationMatrix.h"

float conservidue_pair_score(const Conservidue & c1, const Conservidue & c2) {
	// Simple sum of pairs formula
	float answer = 0.0;

	// nested loops of map iterators seem to cause problems on the Mac
	// So store the c2 keys in a vector
	vector<char> c2_chars;
	map<char,float>::const_iterator res2; // Is the nested loop a problem?
	for (res2 = c2.residue_counts.begin(); 
		 res2 != c2.residue_counts.end();
		 res2 ++) {
		c2_chars.push_back(res2->first);
	}
	
	// nested loops of map iterators seem to cause problems on the Mac
	// So store the c2 keys in a vector
	vector<char> c1_chars;
	map<char,float>::const_iterator res1; // Is the nested loop a problem?
	for (res1 = c1.residue_counts.begin(); 
		 res1 != c1.residue_counts.end();
		 res1 ++) {
		c1_chars.push_back(res1->first);
	}
	
	vector<char>::const_iterator charit1;
	for (charit1 = c1_chars.begin(); charit1 != c1_chars.end(); charit1++) {
			
		vector<char>::const_iterator charit2;
		for (charit2 = c2_chars.begin(); charit2 != c2_chars.end(); charit2++) {
			
			float weight1 = c1.residue_counts.find(*charit1)->second;
			float weight2 = c2.residue_counts.find(*charit2)->second;
		
			float sequence_count = weight1 * weight2;
			float bit_score = blosum62.get_score(*charit1, *charit2);

			answer += sequence_count * bit_score;
		}
	}
	
	return answer;
}

