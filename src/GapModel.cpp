// $Id$
// $Header$
// $Log$
// Revision 1.2  2004/05/28 19:37:31  cmbruns
// Match and extension components of sum-of-pairs are correct in profile alignment.
// Open, close and Delete components are still lower bounds
// New total_penalty subroutine to support sum_of_pairs subroutine
//
// Revision 1.1  2004/05/23 23:46:52  cmbruns
// New structures for piecewise linear gap penalty model
//

#include "GapModel.h"

const GapModel protein_gap_model(PROTEIN_GAP_MODEL_INSTANCE);
int DEFAULT_GAP_SEGMENT_COUNT = protein_gap_model.segment_count();

AlignmentScore GapModel::total_penalty(unsigned int gap_length) const {
	AlignmentScore answer;
	double previous_extension_penalty = 0.0;
	int assigned_gap_length = 0;
	int delta_gap_length = 0;
	for (unsigned int gseg = 0; gseg < segment_count(); gseg++) {
		double current_extension_penalty = extension_penalty(gseg);

		delta_gap_length = minimum_gap_length(gseg) - assigned_gap_length - 1;

		if (minimum_gap_length(gseg) > gap_length) {
			delta_gap_length = gap_length - assigned_gap_length;
			if (delta_gap_length > 0) answer.extension() += delta_gap_length * previous_extension_penalty;
			assigned_gap_length += delta_gap_length;
			break;
		}

		if (delta_gap_length > 0) {
			answer.extension() += delta_gap_length * previous_extension_penalty;
			assigned_gap_length += delta_gap_length;
		}
		
		previous_extension_penalty = current_extension_penalty;
		if (assigned_gap_length < 0) assigned_gap_length = 0;
	}

	delta_gap_length = gap_length - assigned_gap_length;
	if (delta_gap_length > 0) {
		answer.extension() += delta_gap_length * previous_extension_penalty;
		assigned_gap_length += delta_gap_length;
	}

	return answer;
}
