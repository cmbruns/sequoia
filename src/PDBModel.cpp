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
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/05/22 02:14:20  bruns
// Changed FastaSequence to BioSequence, to match new way
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "PDBModel.h"

// PDBModel contains one or more PDBChains
// (for convenience, it also contains its own set of PDBAtoms)
// parent of PDBChain, child of PDBEntry

void PDBModel::check_residues() {
  vector<char>::const_iterator chain_id;
  for (chain_id = chain_order().begin(); 
       chain_id != chain_order().end(); 
       chain_id ++) {
    PDBChain & pdb_chain = chain(*chain_id);
    pdb_chain.check_residues();
  }
}

const vector<char> & PDBModel::chain_order() const {return private_chain_order;}
const map<char, PDBChain> & PDBModel::chain() const {return private_chain;}
void PDBModel::clear() {
  private_chain.clear();
  private_chain_order.clear();
  private_endmdl_card = "";
  private_model_card = "";
}
void PDBModel::set_endmdl_card(string s) {
  private_endmdl_card = s;
}
void PDBModel::set_model_card(string s) {
  private_model_card = s;
}
const string & PDBModel::model_card() const {return private_model_card;}
const string & PDBModel::endmdl_card() const {return private_endmdl_card;}
void PDBModel::add_line(const string & line) {
  if (line.substr(0,6) == "ENDMDL") {
    set_endmdl_card(line);
    return;
  }
  if (line.substr(0,6) == "MODEL ") {
    set_model_card(line);
    return;
  }
  
  // I know how to parse atoms
  if (line.substr(0,6) == "ATOM  ") {
    PDBAtom atom(line);
    PDBChain & pdb_chain = chain(atom.get_chain_id());
    pdb_chain.add_atom(atom);
    return;
  }
  if (line.substr(0,6) == "HETATM") {
    PDBAtom atom(line);
    PDBChain & pdb_chain = chain(atom.get_chain_id());
    pdb_chain.add_atom(atom);
    return;
  }
  if (line.substr(0,6) == "TER   ") {
    PDBAtom fake_atom(line);
    PDBChain & pdb_chain = chain(fake_atom.get_chain_id());
    pdb_chain.add_terminus(fake_atom.get_residue_name(),
		       fake_atom.get_residue_number(),
		       fake_atom.get_insertion_code());
    return;
  }
}

PDBChain & PDBModel::chain(const char chain_id) {
  if (private_chain.find(chain_id) == private_chain.end()) {
    // Make a new chain
    PDBChain new_chain(chain_id);
	  new_chain.title() = get_title();
	  new_chain.sequence_id() = get_pdb_id() + '_' + chain_id;
    private_chain[chain_id] = new_chain;
    private_chain_order.push_back(chain_id);
  }
  return (*(private_chain.find(chain_id))).second;
}

const PDBChain & PDBModel::get_chain(const char chain_id) const {
  return (*(private_chain.find(chain_id))).second;
}

vector<BioSequence> PDBModel::to_fasta(const string & id, const string & title) const {
  vector<BioSequence> fasta;
  // Loop over chains
  vector<char>::const_iterator chain_id;
  for (chain_id = chain_order().begin(); 
       chain_id != chain_order().end(); 
       chain_id ++) {
    const PDBChain & chain = get_chain(*chain_id);
    fasta.push_back(chain.to_fasta(id + "_" + chain.get_chain_id(), title));
  }
  return fasta;
}

vector<BioSequence> PDBModel::to_coord_fasta(const string & id, const string & title) const {
  vector<BioSequence> fasta;
  // Loop over chains
  vector<char>::const_iterator chain_id;
  for (chain_id = chain_order().begin(); 
       chain_id != chain_order().end(); 
       chain_id ++) {
    const PDBChain & chain = get_chain(*chain_id);
    fasta.push_back(chain.to_coord_fasta(id + "_" + chain.get_chain_id(), title));
  }
  return fasta;
}

vector<BioSequence> PDBModel::to_seqres_fasta(const string & id, const string & title) const {
  vector<BioSequence> fasta;
  // Loop over chains
  vector<char>::const_iterator chain_id;
  for (chain_id = chain_order().begin(); 
       chain_id != chain_order().end(); 
       chain_id ++) {
    const PDBChain & chain = get_chain(*chain_id);
    fasta.push_back(chain.to_seqres_fasta(id + "_" + chain.get_chain_id(), title));
  }
  return fasta;
}

ostream & operator<<(ostream & os, const PDBModel & model) {
  if (model.model_card().length() > 0)
    os << model.model_card() << endl;

  vector<char>::const_iterator chain_id;
  for (chain_id = model.chain_order().begin(); chain_id != model.chain_order().end(); ++ chain_id) {
    const PDBChain & chain = model.get_chain(*chain_id);
    os << chain;
  }

  if (model.endmdl_card().length() > 0)
    os << model.endmdl_card() << endl;
  
  return os;
}
