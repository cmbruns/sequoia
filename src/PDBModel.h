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
