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
// Revision 1.4  2004/06/14 16:50:34  cmbruns
// Renamed Score.h to AlignmentScore.h
//
// Revision 1.3  2004/06/04 19:09:06  cmbruns
// Updated GPL header
//
// Change exception syntax to use classes instead of integers
//
// Revision 1.2  2004/05/28 19:39:36  cmbruns
// Match and extension components of sum-of-pairs are correct in profile alignment.
// Open, close and Delete components are still lower bounds
// Changed default protein gap model to have one segment to keep things simple for the moment.
// Temporarily replaced all double alignment scores with AlginmentScore objects for debugging.  NOTE - change these back later for speed and space.
//
// Revision 1.1  2004/05/23 23:46:52  cmbruns
// New structures for piecewise linear gap penalty model
//

// classes and instances for piecewise linear gap penalties

#ifndef __GAP_MODEL_H__
#define __GAP_MODEL_H__

#include <vector>
#include "Exceptions.h"
#include "AlignmentScore.h"

using namespace std;

enum GapModelInstance {
	PROTEIN_GAP_MODEL_INSTANCE
};

// GapSegment is one piecewise component of a GapModel
class GapSegment {
protected:
	unsigned int p_minimum_gap_length;
	AlignmentScore p_extension_penalty;
public:
	GapSegment(unsigned int min_len, double ext_pen) {
		p_minimum_gap_length = min_len;
		p_extension_penalty.extension() = ext_pen;
	}
	unsigned int minimum_gap_length() const {return p_minimum_gap_length;}
	AlignmentScore extension_penalty() const {return p_extension_penalty;}
};

// GapModel encodes a piecewise linear gap extension model
class GapModel {
protected:
	vector<GapSegment> p_segments; // the segments themselves
public:
	unsigned int segment_count() const {return p_segments.size();}
	AlignmentScore extension_penalty(unsigned int segment) const {return p_segments[segment].extension_penalty();}
	AlignmentScore total_penalty(unsigned int gap_length) const;
	unsigned int minimum_gap_length(unsigned int segment) const {return p_segments[segment].minimum_gap_length();}
	void add_segment(unsigned int min_length, double ext_pen) {
		if (ext_pen > 0) throw GAP_MODEL_PROBLEM_EXCEPTION();
		if (p_segments.size() > 0) {
			// Sanity check
			unsigned int previous_segment_min_length = p_segments.back().minimum_gap_length();
			double previous_extension_penalty = p_segments.back().extension_penalty();
			if (ext_pen < previous_extension_penalty) throw GAP_MODEL_PROBLEM_EXCEPTION();
			if (min_length < previous_segment_min_length) throw GAP_MODEL_PROBLEM_EXCEPTION();
		}
		GapSegment new_segment(min_length, ext_pen);
		p_segments.push_back(new_segment);
	}
	
	GapModel(GapModelInstance instance) {
		if (instance == PROTEIN_GAP_MODEL_INSTANCE) {

			// These are derived from FSSP data
			// They seem to be independent of evolutionary distance
			// add_segment(0, -0.90); 
			// add_segment(5, -0.26);
			// add_segment(21, -0.11);
			// add_segment(61, -0.040);
			// add_segment(301, -0.005); // This one is made up

			// Simple model for sum of pairs debugging
			add_segment(0, -1.00);
		}
	}
};

// Gap parameters associated with one residue
class ResidueGapParameter {
	// friend class Conservidue;
protected:
	// Offsets are needed to make later segments worse than earlier ones until their minimum distance comes up
	vector<AlignmentScore> p_open_offset; // initial gap offset for segment
	vector<AlignmentScore> p_close_offset; // final gap offset for segment
	const GapModel * p_gap_model;
	double sequence_weight;
public:
	AlignmentScore open_offset(unsigned int gseg) const {return p_open_offset[gseg];}
	AlignmentScore close_offset(unsigned int gseg) const {return p_close_offset[gseg];}
	AlignmentScore extension_penalty(unsigned int gseg) const {return sequence_weight * p_gap_model->extension_penalty(gseg);}
	// Merge together two profiles during alignment
	ResidueGapParameter combine(const ResidueGapParameter & rgp2) const {
		if (p_gap_model != rgp2.p_gap_model) throw GAP_MODEL_PROBLEM_EXCEPTION();
		ResidueGapParameter answer = *this;
		answer.sequence_weight += rgp2.sequence_weight;
		for (unsigned int gseg = 0; gseg < p_gap_model->segment_count(); gseg ++) {
			answer.p_open_offset[gseg] += rgp2.p_open_offset[gseg];
			answer.p_close_offset[gseg] += rgp2.p_close_offset[gseg];
		}
		return answer;
	}
	unsigned int segment_count() const {return p_gap_model->segment_count();}
	void set_weight(double w) {sequence_weight = w;}
	
	ResidueGapParameter(const GapModel & gap_model, double seq_weight = 0.0) {
		sequence_weight = seq_weight;
		
		p_gap_model = & gap_model;
		// Initialize initial values of open and close offset
		AlignmentScore no_score;
		p_open_offset.assign(p_gap_model->segment_count(), no_score);
		p_close_offset.assign(p_gap_model->segment_count(), no_score);
		unsigned int previous_start = 1;
		double previous_extension_penalty = 0.0;
		double previous_total_penalty = 0.0;
		for (unsigned int gseg = 0; gseg < p_gap_model->segment_count(); gseg ++) {
			unsigned int start = p_gap_model->minimum_gap_length(gseg);
			if (start < 1) start = 1;
			double penalty = p_gap_model->extension_penalty(gseg) * sequence_weight;
			
			// Compute difference between gap penalty up to position start-1, and
			// what it would be if new extension penalty had always been in effect
			previous_total_penalty += (start - previous_start) * previous_extension_penalty;
			double penalty_difference = previous_total_penalty - ((start - 1) * penalty);
			p_open_offset[gseg].extension() = penalty_difference / 2.0;
			p_close_offset[gseg].extension() = penalty_difference / 2.0;
			
			previous_extension_penalty = penalty;
			previous_start = start;
		}
	}
};

extern const GapModel protein_gap_model;
extern int DEFAULT_GAP_SEGMENT_COUNT;

#endif
