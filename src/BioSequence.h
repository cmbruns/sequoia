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
// Revision 1.1  2004/05/11 20:26:12  cmbruns
// Initial revision
//
// Revision 1.8  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.7  2002/09/13 22:51:40  bruns
// Added include <vector> statement
//
// Revision 1.6  2002/05/22 02:06:38  bruns
// Created BioSequence.h by combining FastaSequence.h and BaseSequence.h
//
// Revision 1.5  2002/05/15 19:29:40  bruns
// Replaced operator<< with print subroutine
//
// Revision 1.4  2001/12/18 01:19:19  bruns
// Added new_clone function, required by parent class
//
// Revision 1.3  2001/11/16 00:15:25  bruns
// Broad changes to make it decend from BaseSequence
//
// Revision 1.2  2001/11/15 20:19:55  bruns
// Added cvs tags up through FastaSequence.h
//
#ifndef _BIOSEQUENCE_H_
#define _BIOSEQUENCE_H_

#include "Residue.h"
#include <string>
#include <iostream>
#include <vector>

class SequenceAlignment;

using namespace std;

class BioSequence {
protected:
  vector<Residue *> private_sequence;
  string private_id;
  string private_title;
  BioSequence * new_clone() const {
	  return new BioSequence(*this);
  }
public:

  // default constructor
  BioSequence() {}
  
  // copy constructor
  BioSequence(const BioSequence & s2) {
	  *this = s2;
  }
  
  // simple constructor
  BioSequence(const char * seq_string) {

	  const char * pos = seq_string;
	  int residue_number = 0;
	  while (*pos != NULL) {

		  char one_letter_code = *pos;
		  Residue * residue_pointer = new_protein_residue(one_letter_code);
		  residue_pointer->residue_number() = residue_number;
		  residue_pointer->set_sequence_pointer(this);
		  private_sequence.push_back(residue_pointer);
		  
		  pos ++;
		  residue_number ++;
	  }
  }
  
  // destructor
  ~BioSequence() {clear();}
  
  Macromolecule macromolecule; // TODO -- initialize this
  const string & get_id() const;
  const string & get_title() const;
  const string get_sequence() const {return this->get_string();}
  void set_id(const string & s);
  void set_title(const string & s);
  // void set_sequence(const string & s);

  ostream & print(ostream & os) const;

  // Assignment operator required (just like copy constructor) since 
	// New pointers must be generated for new residues
	BioSequence & operator= (const BioSequence & s2) {
		if (this == &s2) return *this;   // Gracefully handle self assignment
		
		clear(); // delete previous content
		int i;
		for (i = 0; i < s2.length(); ++i) {
			Residue * residue_pointer = s2[i].new_clone();
			residue_pointer->set_sequence_pointer(this);
			private_sequence.push_back(residue_pointer);
		}
		return *this;
	} 

  const string get_string() const;
  // void add_residue(Residue * r) {private_sequence.push_back(r.new_clone());}
  // void add_residue(const char & c) {
  //   private_sequence.push_back(new FastaResidue(c.new_clone()));
  // }
  void add_residue(const Residue & r) {
    private_sequence.push_back(r.new_clone());
  }
  int length() const {return private_sequence.size();}
  const Residue & operator[](int index) const {return *(private_sequence[index]);}
  void clear() {
    vector<Residue *>::iterator i;
    for (i = private_sequence.begin(); i != private_sequence.end(); ++i) {
      delete *i;
    }
    private_sequence.clear();
  }

  bool is_identical_to(const BioSequence &seq2) const;
  bool is_substring_of(const BioSequence &seq2) const;
  bool is_subsequence_of(const BioSequence &seq2) const;
  bool is_erroneous_subsequence_of(const BioSequence &seq2, int max_errors) const;

  bool empty() const;

};

ostream & operator<<(ostream & os, const BioSequence & s);
istream & operator>> (istream & is, BioSequence & s);
// ostream & operator<< (ostream & os, const BioSequence & s);

#endif

