// $Id$
// This file is part of the Sequoia package for macromolecular 
//  sequence/structure analysis
// Copyright (C) 2002  Christopher M. Bruns, Ph.D.
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
// To contact the author, write to cmbruns@attbi.com or bruns@scripps.edu
// In publications please cite: Bruns et al (1999), J.Mol.Biol. 288:427-439
// Please submit bug reports at http://bruns.homeip.net/bugzilla/
// 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.1  2004/05/28 19:32:44  cmbruns
// Match and extension components of sum-of-pairs are correct in profile alignment.
// Open, close and Delete components are still lower bounds
// Changed name from cc to cpp to match other files
//
// Revision 1.3  2004/05/22 19:59:21  cmbruns
// Change name of print routine to print_fasta
// In BioSequence istream<<:
//   change free() to delete() !!!!
//   note explicit gap residues in assigning sequence numbering
//
// Revision 1.2  2004/05/19 00:36:48  cmbruns
// Added print_debug routine for BioSequence object
// Call residue operator<< for output, rather than one_letter_code() routine directly.
//
// Revision 1.1.1.1  2004/05/11 20:26:12  cmbruns
// Initial Max repository for latest sequoia
//
// Revision 1.7  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.6  2002/05/22 02:07:22  bruns
// Created BioSequence.cc by combining FastaSequence.cc with BaseSequence.cc
//
// Revision 1.5  2001/11/19 20:44:46  bruns
// Changed BioSequence::add_residue calls to take character argument
//
// Revision 1.4  2001/11/16 01:07:59  bruns
// Fixed bug 5, so that BioSequence stores an array of pointers to Residue
// objects allocated with new.
//
// Revision 1.3  2001/11/16 00:14:32  bruns
// Broad changes to make it compile as descendant of BioSequence
//
// Revision 1.2  2001/11/15 20:19:55  bruns
// Added cvs tags up through BioSequence.h
//
#include "BioSequence.h"
#include "Residue.h"
#include "SequenceAlignment.h"

#define MAX_LINE_LENGTH 10000

int CHARACTERS_PER_LINE = 60; // in fasta output

const Residue & BioSequence::operator[](unsigned int i) const {return *private_sequence[i];}
const string BioSequence::get_sequence() const {return this->get_string();}
void BioSequence::add_residue(const Residue & r) { // residue allocated with new
	private_sequence.push_back(r.new_clone());
	Residue & residue = *private_sequence.back();
	residue.set_sequence_pointer(this);
	residue.sequence_residue_index = private_sequence.size() - 1;
}
unsigned int BioSequence::length() const {return private_sequence.size();}
void BioSequence::clear() {
    vector<Residue *>::iterator i;
    for (i = private_sequence.begin(); i != private_sequence.end(); ++i) {
		delete *i;
    }
    private_sequence.clear();
}

const string & BioSequence::get_id() const {return private_id;}
const string & BioSequence::get_title() const {return private_title;}
void BioSequence::set_id(const string & s) {private_id = s;}
void BioSequence::set_title(const string & s) {private_title = s;}

ostream & BioSequence::print_debug(ostream & os, unsigned int indent_size) const {
	string indent = "";
	for(unsigned int i=0;i<indent_size;i++)indent += " ";
	
	os << indent << "sequence pointer = " << this << endl;
	os << indent << "Residues:" << endl;
	for (unsigned int i = 0; i < private_sequence.size(); i++) {
		private_sequence[i]->print_debug(os, indent_size + 2);
	}
	os << endl;

	return os;
}

ostream & BioSequence::print_fasta(ostream & os) const {
    const BioSequence & sequence = *this;

    os << ">";
    os << get_id();
    os << "\t";
    os << get_title();
    os << endl;
    for (unsigned int i = 0; i < length(); ++i) {
		if ((i > 0) && (i % CHARACTERS_PER_LINE == 0))
			os << endl; // end of line every CHARACTERS_PER_LINE characters
		os << sequence[i];
    }
    os << endl;
    return os;
}

istream & BioSequence::read_fasta(istream & is) {
	BioSequence & s = *this;

	char this_char = '\n';
	char line_buffer[MAX_LINE_LENGTH + 1];

  s.set_id("");
  s.set_title("");
  s.clear();

  // Get one line at a time, checking for '>' at beginning
  // 1 look for '>'
  while ((this_char != '>') && is.good()) {
    if (this_char != '\n') // get rest of line if not at start of line
      is.getline(line_buffer, MAX_LINE_LENGTH);
    this_char = is.get();
  }
  if (!is.good()) return is;
  is.getline(line_buffer, MAX_LINE_LENGTH);
  string header = line_buffer;
  unsigned int id_start = header.find_first_not_of(" \t\n\v\b\r\f");
  if (id_start == string::npos) {
    s.set_id("");
    s.set_title("");
  } else {
    unsigned int id_end = header.find_first_of(" \t\n\v\b\r\f", id_start);
    if (id_end == string::npos) {
      s.set_id(header.substr(id_start));
      s.set_title("");
    } else {
      s.set_id(header.substr(id_start, id_end - id_start));
      unsigned int title_start = header.find_first_not_of(" \t\n\v\b\r\f", id_end);
      if (title_start == string::npos) 
	s.set_title("");
      else
	s.set_title(header.substr(title_start));
    }
  }
  
  // 3 eat sequence until next '\n>'
  int residue_number = 0;
  s.clear();
  this_char = is.get();
  while ((this_char != '>') && is.good()) {
    if (this_char != '\n') { // get rest of line if not at start of line

		Residue * res_ptr = new_protein_residue(this_char);
		if (res_ptr->is_gap()) // Gap
			res_ptr->residue_number() = -1;
		else { // Not a gap
			residue_number ++;
			res_ptr->residue_number() = residue_number;
		}
		s.add_residue(*res_ptr);
		delete res_ptr;

		is.getline(line_buffer, MAX_LINE_LENGTH);
		int i;
		for (i = 0; line_buffer[i] != '\0'; ++i) {

			if (i >= MAX_LINE_LENGTH) break;

			res_ptr = new_protein_residue(line_buffer[i]);
			if (res_ptr->is_gap()) // Gap
				res_ptr->residue_number() = -1;
			else { // Not a gap
				residue_number ++;
				res_ptr->residue_number() = residue_number;
			}
			s.add_residue(*res_ptr);
			delete res_ptr;
		}
    }
    this_char = is.get();
  }
  // Put '>' back on the stream
  if (this_char == '>') is.putback(this_char);
  return is;
}

bool BioSequence::is_identical_to(const BioSequence &seq2) const {
  return (get_string() == seq2.get_string());
}

bool BioSequence::is_substring_of(const BioSequence &seq2) const {
  if (length() > seq2.length()) return false;
  unsigned int pos = seq2.get_string().find(get_string());
  if (pos == string::npos) return false;
  return true;
}

bool BioSequence::is_subsequence_of(const BioSequence &seq2) const {
  int l1 = length();
  int l2 = seq2.length();
  if (l1 > l2) return false;

  const BioSequence & seq1 = *this;
  int pos1, pos2;
  for (pos1 = 0, pos2 = 0; 
       pos1 < l1; 
       pos1 ++, pos2 ++) {
    char r1 = (seq1[pos1]).one_letter_code();
    char r2 = (seq2[pos2]).one_letter_code();
    while (r1 != r2) {
      pos2 ++;
      // Are there enough characters left in string 2?
      if ((l2 - pos2) < (l1 - pos1)) return false;
      r2 = (seq2[pos2]).one_letter_code();
    }
  }

  return true;
}

const string BioSequence::get_string() const {
	string answer = "";
	vector<Residue *>::const_iterator i;
	const vector<Residue *> & v = private_sequence;
	for (i = v.begin(); i != v.end(); ++i) {
		answer += (*i)->one_letter_code();
	}
	return answer;
}

bool BioSequence::is_erroneous_subsequence_of(const BioSequence &seq2, int max_errors) const {
  const BioSequence & seq1 = *this;

  int l1 = length();
  int l2 = seq2.length();
  if (l1 > l2) return false;

  int pos2[max_errors + 1]; // earliest position in seq2 with given number of errors

  // Simulate a fancy "for" loop
  // Initialization condition
  int pos1 = 0; // position in sequence 1
  int e;
  for (e = 0; e <= max_errors; ++e) pos2[e] = 0;
  while (1) {
    char r1 = (seq1[pos1]).one_letter_code();
    for (e = max_errors; e >= 0; --e) {
      if (pos2[e] >= l2) continue;
      char r2 = (seq2[pos2[e]]).one_letter_code();

      // Compare 3 cases
      //  1 - no new errors; march from previous errors = e position
      //  2 - new error at current position
      //  3 - insertion at current position - FIXME

      // Case 1 - continue without new errors
      while (r1 != r2) {
	pos2[e] ++;
	// Are there enough characters left in string 2?
	if (pos2[e] >= l2) break;
	if ((l2 - pos2[e]) < (l1 - pos1)) break;
	r2 = (seq2[pos2[e]]).one_letter_code();
      }

      // Case 2 - new error position, extend from e-1 error position
      if ((e > 0) && (pos2[e - 1] < pos2[e])) {
	pos2[e] = pos2[e - 1];
      }

    }
    // Increment condition
    pos1 ++;
    for(e = 0; e <= max_errors; ++e)
      pos2[e] ++;
    // Termination condition
    if (pos1 >= l1) break;
  }
  if (pos2[max_errors] <= l2) return true;
  else return false;
}

bool BioSequence::is_empty() const {
  return get_string().empty();
}

// Assignment operator required (just like copy constructor) since 
// New pointers must be generated for new residues
BioSequence & BioSequence::operator= (const BioSequence & s2) {
	if (this == &s2) return *this;   // Gracefully handle self assignment
	
	clear(); // delete previous content
	for (unsigned int i = 0; i < s2.length(); ++i) {
		add_residue(s2[i]);
	}
	
	private_id = s2.private_id;
	private_title = s2.private_title;
	macromolecule = s2.macromolecule;
	alignment_sequence_index = s2.alignment_sequence_index;
	parent_alignment = NULL;
	p_gap_model = s2.p_gap_model;
	left_gap_factor = s2.left_gap_factor;
	right_gap_factor = s2.right_gap_factor;
	left_gap_extension_factor = s2.left_gap_extension_factor;
	right_gap_extension_factor = s2.right_gap_extension_factor;
	set_weight(s2.get_weight());
	
	return *this;
} 

// default constructor
BioSequence::BioSequence() : 
  private_id("(no ID)"), 
  private_title("(no description)"),
  macromolecule(PROTEIN_MOLECULE),
  alignment_sequence_index(-1),
  parent_alignment(NULL),
  p_gap_model(& protein_gap_model),
  left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
  right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
  left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
  right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
  p_weight(1.0)
{}

// copy constructor
BioSequence::BioSequence(const BioSequence & s2) {
	*this = s2;
}

// BioSequence string constructor
BioSequence::BioSequence(const char * seq_string)  : 
  private_id("(no ID)"), 
  private_title("(no description)"),
  macromolecule(PROTEIN_MOLECULE),
  alignment_sequence_index(-1),
  parent_alignment(NULL),
  p_gap_model(& protein_gap_model),
  left_gap_factor(DEFAULT_LEFT_GAP_FACTOR),
  right_gap_factor(DEFAULT_RIGHT_GAP_FACTOR),
  left_gap_extension_factor(DEFAULT_LEFT_GAP_EXTENSION_FACTOR),
  right_gap_extension_factor(DEFAULT_RIGHT_GAP_EXTENSION_FACTOR),
  p_weight(1.0)
{
	
	const char * pos = seq_string;
	int residue_number = 0;
	while (*pos != NULL) {
		
		char one_letter_code = *pos;
		Residue * residue_pointer = new_protein_residue(one_letter_code);
		if (!residue_pointer->is_gap()) {
			residue_number ++;
			residue_pointer->residue_number() = residue_number;
		}
		add_residue(*residue_pointer);
		delete residue_pointer;
		
		pos ++;
	}
}

// destructor
BioSequence::~BioSequence() {
	clear();
}

istream & operator>> (istream & is, BioSequence & s) {
	return s.read_fasta(is);
}

ostream & operator<<(ostream & os, const BioSequence & s) {
  s.print_fasta(os);
  return os;
}

