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

// $Id$
// $Header$
// $Log$
// Revision 1.2  2004/06/14 16:38:00  cmbruns
// Added minus operator
// renamed AlignmentScore to DebugAlignmentScore
// sketched mechanism for reverting AlignmentScore to double once testing is complete
//
// Revision 1.1  2004/06/07 18:47:02  cmbruns
// Renamed Score object to AlignmentScore
//
// Revision 1.4  2004/06/04 19:15:46  cmbruns
// Updated GPL header
//
// Moved conservidue pair score function from here to be a member function of conservidue
//

#ifndef __ALIGNMENT_SCORE_H__
#define __ALIGNMENT_SCORE_H__

#include <iostream>

using namespace std;

// Slow score object for easire debugging
#define __DEBUG_ALIGNMENT_SCORES__ 1
#ifdef __DEBUG_ALIGNMENT_SCORES__
#define AlignmentScore DebugAlignmentScore
#else
#define AlignmentScore double
#define .match() ""
#define .extension() ""
#define .opening() ""
#define .closing() ""
#define .deletion() ""
#define .transition() ""
#define ->match() ""
#define ->extension() ""
#define ->opening() ""
#define ->closing() ""
#define ->deletion() ""
#define ->transition() ""
#endif

enum DebugAlignmentScoreType {
	MATCH_SCORE,
	EXTENSION_SCORE,
	OPENING_SCORE,
	CLOSING_SCORE,
	DELETION_SCORE,
	TRANSITION_SCORE
};

// Capture details of alignment score
// Want to be able to replace this class with double for efficiency
class DebugAlignmentScore {
	
private:
	
	double p_match_score; // contribution from aligned scores
	double p_extension_score; // contribution from aligned scores
	double p_opening_score; // contribution from aligned scores
	double p_closing_score; // contribution from aligned scores
	double p_deletion_score; // contribution from aligned scores
	double p_transition_score; // contribution from aligned scores
	
public:
	operator double() const {return get_total();}
	void clear() {
		match() = 0;
		extension() = 0;
		opening() = 0;
		closing() = 0;
		deletion() = 0;
		transition() = 0;
	}
	
	double get_match() const {return p_match_score;}
	double get_extension() const {return p_extension_score;}
	double get_opening() const {return p_opening_score;}
	double get_closing() const {return p_closing_score;}
	double get_deletion() const {return p_deletion_score;}
	double get_transition() const {return p_transition_score;}
	
	double & match() {return p_match_score;}
	double & extension() {return p_extension_score;}
	double & opening() {return p_opening_score;}
	double & closing() {return p_closing_score;}
	double & deletion() {return p_deletion_score;}
	double & transition() {return p_transition_score;}
	
	double get_total() const {return 
		get_match() + 
		get_extension() + 
		get_opening() +
		get_closing() +
		get_deletion() +
		get_transition();}
	DebugAlignmentScore operator+(const DebugAlignmentScore & score2) const {
		DebugAlignmentScore answer = *this;
		answer.match() += score2.get_match();
		answer.extension() += score2.get_extension();
		answer.opening() += score2.get_opening();
		answer.closing() += score2.get_closing();
		answer.deletion() += score2.get_deletion();
		answer.transition() += score2.get_transition();
		return answer;
	}
	DebugAlignmentScore operator-(const DebugAlignmentScore & score2) const {
		DebugAlignmentScore answer = *this;
		answer.match() -= score2.get_match();
		answer.extension() -= score2.get_extension();
		answer.opening() -= score2.get_opening();
		answer.closing() -= score2.get_closing();
		answer.deletion() -= score2.get_deletion();
		answer.transition() -= score2.get_transition();
		return answer;
	}
	DebugAlignmentScore & operator+=(const DebugAlignmentScore & score2)  {
		*this = *this + score2;
		return *this;
	}
	DebugAlignmentScore & operator*=(const double r)  {
		*this = *this * r;
		return *this;
	}
	DebugAlignmentScore operator*(const double r) const {
		DebugAlignmentScore answer = *this;
		answer.match() *= r;
		answer.extension() *= r;
		answer.opening() *= r;
		answer.closing() *= r;
		answer.deletion() *= r;
		answer.transition() *= r;
		return answer;
	}
	ostream & print_details(ostream & os = cout) const {
		os << "match score = " << get_match() << endl;
		os << "extension score = " << get_extension() << endl;
		os << "opening score = " << get_opening() << endl;
		os << "closing score = " << get_closing() << endl;
		os << "deletion score = " << get_deletion() << endl;
		os << "transition score = " << get_transition() << endl;
		os << "total score = " << get_total() << endl;
		return os;
	}
	
	DebugAlignmentScore() {
		clear();
	}
	DebugAlignmentScore(DebugAlignmentScoreType score_type, double value) {
		clear();
		switch (score_type) {
			case MATCH_SCORE:
				match() = value;
				break;
			case EXTENSION_SCORE:
				extension() = value;
				break;
			case OPENING_SCORE:
				opening() = value;
				break;
			case CLOSING_SCORE:
				closing() = value;
				break;
			case DELETION_SCORE:
				deletion() = value;
				break;
			case TRANSITION_SCORE:
				transition() = value;
				break;
			default:
				break;
		}
	}	
};

DebugAlignmentScore operator*(const double r, const DebugAlignmentScore & score);

class Conservidue;

#endif
