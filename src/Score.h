#ifndef __SCORE_H__
#define __SCORE_H__

#include <iostream>

using namespace std;

enum AlignmentScoreType {
	MATCH_SCORE,
	EXTENSION_SCORE,
	OPENING_SCORE,
	CLOSING_SCORE,
	DELETION_SCORE,
	TRANSITION_SCORE
};

// Capture details of alignment score
// Want to be able to replace this class with double for efficiency
class AlignmentScore {
	
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
	AlignmentScore operator+(const AlignmentScore & score2) const {
		AlignmentScore answer = *this;
		answer.match() += score2.get_match();
		answer.extension() += score2.get_extension();
		answer.opening() += score2.get_opening();
		answer.closing() += score2.get_closing();
		answer.deletion() += score2.get_deletion();
		answer.transition() += score2.get_transition();
		return answer;
	}
	AlignmentScore & operator+=(const AlignmentScore & score2)  {
		*this = *this + score2;
		return *this;
	}
	AlignmentScore & operator*=(const double r)  {
		*this = *this * r;
		return *this;
	}
	AlignmentScore operator*(const double r) const {
		AlignmentScore answer = *this;
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
	
	AlignmentScore() {
		clear();
	}
	AlignmentScore(AlignmentScoreType score_type, double value) {
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

AlignmentScore operator*(const double r, const AlignmentScore & score);

class Conservidue;

// return bit score of aligning two conservidues
AlignmentScore conservidue_pair_score(const Conservidue & c1, const Conservidue & c2);

#endif
