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
// 

// $Id$
// $Header$
// $Log$
// Revision 1.3  2004/06/04 19:08:05  cmbruns
// Updated GPL header
//
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
