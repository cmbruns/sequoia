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
// Revision 1.1  2004/06/04 19:34:47  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.5  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.4  2002/05/22 02:18:19  bruns
// Changed FastaSequence to BioSequence, to match new way
//
// Revision 1.3  2002/05/15 19:25:35  bruns
// Added get_first_chain, get_first_model, and get_first_residue methods
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _PDBENTRY_H_
#define _PDBENTRY_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "PDBModel.h"
#include "Seqres.h"
#include "Modres.h"

// Object corresponding to a single PDB file

// This is the top level object of this Hierarchy:
// PDBEntry - object representing one PDB file
//  PDBModel - one "model" in a PDB file (e.g. among multiple NMR models)
//   PDBChain - one polymer or protein or molecule
//    PDBResidue - one building block (e.g. one amino acid residue)
//     PDBAtom - one atom

class PDBEntry;
istream & operator>>(istream & input_stream, PDBEntry & pdb);
ostream & operator<<(ostream & os, const PDBEntry & pdb);

// Parent of PDBModel
// One PDBEntry contains one or more PDBModels
// Usually only NMR entries have more than one PDBModel
class PDBEntry {
private:
  vector<PDBModel> private_model; // Child coordinate models
  vector<Modres> private_modres; // annotation of modified residues
  map<char,Seqres> private_seqres; // explicit molecular sequence per chain
  string private_id_code;
  string private_classification;
  string private_deposition_date;
  string private_title;
  bool read_header_only;
  bool read_first_model_only;
public:
  PDBEntry();
  void set_header_only(bool flag);
  void set_first_model_only(bool flag);
  bool get_header_only();
  bool get_first_model_only();
  const map<char,Seqres> & get_seqres() const;
  map<char,Seqres> & seqres();
  const string & get_id_code() const;
  const string & get_title() const;
  string & id_code();
  string & classification();
  string & deposition_date();
  string & title();
  vector<Modres> & modres();
  void add_model(PDBModel & model);

  vector<PDBModel> & model();
  const vector<PDBModel> & get_model() const;
  PDBModel & get_first_model() {return model()[0];}
  PDBChain & get_first_chain() {return get_first_model().get_first_chain();}
  PDBResidue & get_first_residue() {return get_first_chain().get_first_residue();}

  vector<BioSequence> to_fasta() const;
  vector<BioSequence> to_seqres_fasta() const;
  vector<BioSequence> to_coord_fasta() const;

  void load_pdb_file(const char * file_name) {
	  ifstream infile(file_name);
	  if (infile == 0) {
		  cerr << "*** ERROR *** : Unable to open file " << file_name << endl;
		  throw NO_SUCH_FILE_EXCEPTION();
	  }
	  infile >> *this;
	  infile.close();
  }
  
  void clear();
};


#endif
