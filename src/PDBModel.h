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
// Revision 1.1  2004/06/04 19:34:48  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.5  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.4  2002/05/22 02:14:04  bruns
// Changed FastaSequence to BioSequence, to match new way
//
// Revision 1.3  2002/05/15 19:24:21  bruns
// Added get_first_chain() method
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _PDBMODEL_H_
#define _PDBMODEL_H_

// Object corresponding to a single model in a PDB file

#include <iostream>
#include <vector>
#include <map>
#include "PDBAtom.h"
#include "PDBChain.h"
#include "BioSequence.h"

// PDBModel contains one or more PDBChains
// (for convenience, it also contains its own set of PDBAtoms)
// parent of PDBChain, child of PDBEntry
class PDBModel {
private:
  map<char, PDBChain> private_chain;
  vector<char> private_chain_order;
  string private_endmdl_card;
  string private_model_card;
  string p_pdb_id;
  string p_title;
public:
  void set_pdb_id(string i) {
	  p_pdb_id = i;
	  map<char, PDBChain>::iterator c;
	  for (c = private_chain.begin(); c != private_chain.end(); c++) {
		  char chain_id = c->first;
		  PDBChain & chain = c->second;
		  chain.sequence_id() = get_pdb_id() + '_' + chain_id;
	  }
  }
  void set_title(string t) {
	  p_title = t;
	  map<char, PDBChain>::iterator c;
	  for (c = private_chain.begin(); c != private_chain.end(); c++) {
		  PDBChain & chain = c->second;
		  chain.title() = get_title();
	  }
  }
  string get_title() {return p_title;}
  string get_pdb_id() {return p_pdb_id;}
  
  const vector<char> & chain_order() const;
  const map<char, PDBChain> & chain() const;
  void clear();

  void set_endmdl_card(string s);
  void set_model_card(string s);

  const string & model_card() const;
  const string & endmdl_card() const;

  const PDBChain & get_chain(const char chain_id) const;
  PDBChain & chain(const char chain_id);
  PDBChain & get_first_chain() {return chain(chain_order()[0]);}

  void add_line(const string & line);

  vector<BioSequence> to_fasta(const string & id, const string & title) const;
  vector<BioSequence> to_coord_fasta(const string & id, const string & title) const;
  vector<BioSequence> to_seqres_fasta(const string & id, const string & title) const;

  void check_residues();
};

ostream & operator<<(ostream & os, const PDBModel & model);

#endif
