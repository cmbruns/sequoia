#ifndef __GAP_MODEL_H__
#define __GAP_MODEL_H__

// classes and instances for piecewise linear gap penalties

// $Id$
// $Header$
// $Log$
// Revision 1.1  2004/05/23 23:46:52  cmbruns
// New structures for piecewise linear gap penalty model
//

using namespace std;
#include <vector>
#include "Exceptions.h"

enum GapModelInstance {
	PROTEIN_GAP_MODEL_INSTANCE
};

// GapSegment is one piecewise component of a GapModel
class GapSegment {
protected:
	unsigned int p_minimum_gap_length;
	double p_extension_penalty;
public:
	GapSegment(unsigned int min_len, double ext_pen) {
		p_minimum_gap_length = min_len;
		p_extension_penalty = ext_pen;
	}
	unsigned int minimum_gap_length() const {return p_minimum_gap_length;}
	double extension_penalty() const {return p_extension_penalty;}
};

// GapModel encodes a piecewise linear gap extension model
class GapModel {
protected:
	vector<GapSegment> p_segments; // the segments themselves
public:
	unsigned int segment_count() const {return p_segments.size();}
	double extension_penalty(unsigned int segment) const {return p_segments[segment].extension_penalty();}
	unsigned int minimum_gap_length(unsigned int segment) const {return p_segments[segment].minimum_gap_length();}
	void add_segment(unsigned int min_length, double ext_pen) {
		if (ext_pen > 0) throw GAP_MODEL_PROBLEM_EXCEPTION;
		if (p_segments.size() > 0) {
			// Sanity check
			unsigned int previous_segment_min_length = p_segments.back().minimum_gap_length();
			double previous_extension_penalty = p_segments.back().extension_penalty();
			if (ext_pen < previous_extension_penalty) throw GAP_MODEL_PROBLEM_EXCEPTION;
			if (min_length < previous_segment_min_length) throw GAP_MODEL_PROBLEM_EXCEPTION;
		}
		GapSegment new_segment(min_length, ext_pen);
		p_segments.push_back(new_segment);
	}
	
	GapModel(GapModelInstance instance) {
		if (instance == PROTEIN_GAP_MODEL_INSTANCE) {
			// These are derived from FSSP data
			// They seem to be independent of evolutionary distance
			add_segment(0, -0.90); 
			add_segment(5, -0.26);
			add_segment(21, -0.11);
			add_segment(61, -0.040);
			add_segment(301, -0.005); // This one is made up
		}
	}
};

// Gap parameters associated with one residue
class ResidueGapParameter {
protected:
	// Offsets are needed to make later segments worse than earlier ones until their minimum distance comes up
	vector<double> p_open_offset; // initial gap offset for segment
	vector<double> p_close_offset; // final gap offset for segment
	const GapModel * p_gap_model;
	double sequence_weight;
public:
	double open_offset(unsigned int gseg) const {return p_open_offset[gseg];}
	double close_offset(unsigned int gseg) const {return p_close_offset[gseg];}
	double extension_penalty(unsigned int gseg) const {return sequence_weight * p_gap_model->extension_penalty(gseg);}
	// Merge together two profiles during alignment
	ResidueGapParameter combine(const ResidueGapParameter & rgp2) const {
		if (p_gap_model != rgp2.p_gap_model) throw GAP_MODEL_PROBLEM_EXCEPTION;
		ResidueGapParameter answer = *this;
		answer.sequence_weight += rgp2.sequence_weight;
		for (unsigned int gseg = 0; gseg < p_gap_model->segment_count(); gseg ++) {
			answer.p_open_offset[gseg] += rgp2.p_open_offset[gseg];
			answer.p_close_offset[gseg] += rgp2.p_close_offset[gseg];
		}
		return answer;
	}
	ResidueGapParameter(const GapModel & gap_model, double seq_weight = 1.0) {
		sequence_weight = seq_weight;
		p_gap_model = & gap_model;
		// Initialize initial values of open and close offset
		p_open_offset.assign(p_gap_model->segment_count(), 0);
		p_close_offset.assign(p_gap_model->segment_count(), 0);
		unsigned int previous_start = 1;
		double previous_extension_penalty = 0.0;
		double previous_total_penalty = 0.0;
		for (unsigned int gseg = 0; gseg < p_gap_model->segment_count(); gseg ++) {
			unsigned int start = p_gap_model->minimum_gap_length(gseg);
			if (start < 1) start = 1;
			double penalty = p_gap_model->extension_penalty(gseg);
			
			// Compute difference between gap penalty up to position start-1, and
			// what it would be if new extension penalty had always been in effect
			previous_total_penalty += (start - previous_start) * previous_extension_penalty;
			double penalty_difference = previous_total_penalty - ((start - 1) * penalty);
			p_open_offset[gseg] = penalty_difference / 2.0;
			p_close_offset[gseg] = penalty_difference / 2.0;
			
			previous_extension_penalty = penalty;
			previous_start = start;
		}
	}
};

extern const GapModel protein_gap_model;

#endif
