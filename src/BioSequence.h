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
// 
// $Id$
//
// $Header$
//
// $Log$
// Revision 1.5  2004/06/04 19:02:35  cmbruns
// Updated GPL header
//
// Revision 1.4  2004/05/28 19:36:29  cmbruns
// Match and extension components of sum-of-pairs are correct in profile alignment.
// Open, close and Delete components are still lower bounds
// Moved most function bodies to .cpp file
// Moved end gap parameters from SequenceAlignment, here to BioSequence
//
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

#include <string>
#include <iostream>
#include <vector>
#include "Residue.h"
#include "GapModel.h"

#define DEFAULT_LEFT_GAP_FACTOR 0.5
#define DEFAULT_RIGHT_GAP_FACTOR 0.5
#define DEFAULT_LEFT_GAP_EXTENSION_FACTOR 0.1
#define DEFAULT_RIGHT_GAP_EXTENSION_FACTOR 0.1

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
  const GapModel * p_gap_model; // gap extension scores

  // Scale scores of end gaps by this amount, but not extensions
  double left_gap_factor; // coefficient for left end gaps, other than extensions
  double right_gap_factor; // coefficient for right end gaps, other than extensions
						   // Additional scale factor applies only to extension portion of gap penalty
  double left_gap_extension_factor;
  double right_gap_extension_factor;	
  
  float p_weight; // Relative weight in sequence alignment

  BioSequence * new_clone() const {
	  return new BioSequence(*this);
  }
public:
  void add_residue(const Residue & r);
  void clear();
  const string & get_id() const;
  const string get_sequence() const;
  const string get_string() const;
  const string & get_title() const;
  float get_weight() const {return p_weight;}
  bool is_empty() const;
  bool is_erroneous_subsequence_of(const BioSequence &seq2, int max_errors) const;
  bool is_identical_to(const BioSequence &seq2) const;
  bool is_substring_of(const BioSequence &seq2) const;
  bool is_subsequence_of(const BioSequence &seq2) const;
  unsigned int length() const;
  const Residue & operator[](unsigned int i) const;
  ostream & print_debug(ostream & os, unsigned int indent_size = 0) const;
  ostream & print_fasta(ostream & os) const;
  istream & read_fasta(istream & is);
  void set_id(const string & s);
  void set_title(const string & s);
  void set_weight(float w) {p_weight = w;}
  // float & weight() {return p_weight;}

  BioSequence & operator= (const BioSequence & s2);
  BioSequence();
  BioSequence(const BioSequence & s2);
  BioSequence(const char * seq_string);
  ~BioSequence();
};

ostream & operator<<(ostream & os, const BioSequence & s);
istream & operator>>(istream & is, BioSequence & s);

#endif

