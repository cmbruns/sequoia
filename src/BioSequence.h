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
// Revision 1.3  2004/05/23 23:45:11  cmbruns
// renamed print to print_fasta
//
// Revision 1.2  2004/05/19 00:55:02  cmbruns
// Closer integration of alignment, sequence, conservidue and residue:
//   added friend classes to BioSequence
//   added parent_alignment pointer
//   added alignment_sequence_index
// moved macromolecule to protected from public
// added initialization block to constructors
// enhanced add_residue subroutine
// used add_residue for all residue creation
// fixed memory leak in residue allocation (delete)
// add operator[]
// add print_debug function
//
// Revision 1.1.1.1  2004/05/11 20:26:12  cmbruns
// Initial Max repository for latest sequoia
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
	friend class SequenceAlignment;
	friend class Conservidue;
protected:
  vector<Residue *> private_sequence; // These are pointers, but we need to delete the sequences
  string private_id;
  string private_title;
  Macromolecule macromolecule; // TODO -- initialize this
  int alignment_sequence_index; // my placement in containing alignment, if any
  SequenceAlignment * parent_alignment;

  BioSequence * new_clone() const {
	  return new BioSequence(*this);
  }
public:

	// default constructor
	BioSequence() : 
		private_id("(no ID)"), 
		private_title("(no description)"), 
		alignment_sequence_index(-1),
		parent_alignment(NULL)
	{}
  
  // copy constructor
  BioSequence(const BioSequence & s2) {
	  *this = s2;
  }
  
  // simple constructor
  BioSequence(const char * seq_string)  : 
		private_id("(no ID)"), 
		private_title("(no description)"), 
		alignment_sequence_index(-1),
		parent_alignment(NULL)
	{

	  const char * pos = seq_string;
	  int residue_number = 0;
	  while (*pos != NULL) {

		  char one_letter_code = *pos;
		  Residue * residue_pointer = new_protein_residue(one_letter_code);
		  add_residue(*residue_pointer);
		  delete residue_pointer;
		  
		  pos ++;
		  residue_number ++;
	  }
  }
  
  // destructor
  ~BioSequence() {
	  clear();
  }
  
  const string & get_id() const;
  const string & get_title() const;
  const Residue & operator[](unsigned int i) const {return *private_sequence[i];}
  const string get_sequence() const {return this->get_string();}
  void set_id(const string & s);
  void set_title(const string & s);
  // void set_sequence(const string & s);

  ostream & print_debug(ostream & os, unsigned int indent_size = 0) const;
  ostream & print_fasta(ostream & os) const;

  // Assignment operator required (just like copy constructor) since 
	// New pointers must be generated for new residues
	BioSequence & operator= (const BioSequence & s2) {
		if (this == &s2) return *this;   // Gracefully handle self assignment
		
		clear(); // delete previous content
		for (unsigned int i = 0; i < s2.length(); ++i) {
			add_residue(s2[i]);
		}
		
		private_id = s2.private_id;
		private_title = s2.private_title;
		macromolecule = s2.macromolecule;
		
		return *this;
	} 

  const string get_string() const;
  // void add_residue(Residue * r) {private_sequence.push_back(r.new_clone());}
  // void add_residue(const char & c) {
  //   private_sequence.push_back(new FastaResidue(c.new_clone()));
  // }
  void add_residue(const Residue & r) { // residue allocated with new
	  private_sequence.push_back(r.new_clone());
	  Residue & residue = *private_sequence.back();
	  residue.set_sequence_pointer(this);
	  residue.sequence_residue_index = private_sequence.size() - 1;
  }
  unsigned int length() const {return private_sequence.size();}
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

#endif

