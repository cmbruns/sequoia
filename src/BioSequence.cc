/* Copyright (c) 2005 Christopher M. Bruns
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

// 
// $Id$
//
// $Header$
//
// $Log$
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

const string & BioSequence::get_id() const {return private_id;}
const string & BioSequence::get_title() const {return private_title;}
void BioSequence::set_id(const string & s) {private_id = s;}
void BioSequence::set_title(const string & s) {private_title = s;}
// void BioSequence::set_sequence(const string & s) {private_sequence = s;}

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

istream & operator>> (istream & is, BioSequence & s) {
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

const string BioSequence::get_string() const {
  string answer = "";
  vector<Residue *>::const_iterator i;
  const vector<Residue *> & v = private_sequence;
  for (i = v.begin(); i != v.end(); ++i) {
    answer += (*i)->one_letter_code();
  }
  return answer;
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

bool BioSequence::empty() const {
  return get_string().empty();
}

ostream & operator<<(ostream & os, const BioSequence & s) {
  s.print_fasta(os);
  return os;
}

