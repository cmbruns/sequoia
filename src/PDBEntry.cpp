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
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/05/22 02:18:19  bruns
// Changed FastaSequence to BioSequence, to match new way
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "PDBEntry.h"

// Object corresponding to a single PDB file

// This is the top level object of this Hierarchy:
// PDBEntry - object representing one PDB file
//  PDBModel - one "model" in a PDB file (e.g. among multiple NMR models)
//   PDBChain - one polymer or protein or molecule
//    PDBResidue - one building block (e.g. one amino acid residue)
//     PDBAtom - one atom

PDBEntry::PDBEntry() {
  set_header_only(false);
  set_first_model_only(false);
}
void PDBEntry::set_header_only(bool flag) {
  read_header_only = flag;
}
void PDBEntry::set_first_model_only(bool flag) {
  read_first_model_only = flag;
}
bool PDBEntry::get_header_only() {
  return read_header_only;
}
bool PDBEntry::get_first_model_only() {
  return read_first_model_only;
}
const map<char,Seqres> & PDBEntry::get_seqres() const {
  return private_seqres;
}
map<char,Seqres> & PDBEntry::seqres() {
  return private_seqres;
}
const string & PDBEntry::get_id_code() const {return private_id_code;}
const string & PDBEntry::get_title() const {return private_title;}
string & PDBEntry::id_code() {return private_id_code;}
string & PDBEntry::classification() {return private_classification;}
string & PDBEntry::deposition_date() {return private_deposition_date;}
string & PDBEntry::title() {return private_title;}
vector<PDBModel> & PDBEntry::model() {return private_model;}
const vector<PDBModel> & PDBEntry::get_model() const {return private_model;}
vector<Modres> & PDBEntry::modres() {return private_modres;}
void PDBEntry::add_model(PDBModel & model) {
	private_model.push_back(model);
	PDBModel & m = private_model.back();
	m.set_title(title());
	m.set_pdb_id(id_code());
}

vector<BioSequence> PDBEntry::to_fasta() const {
  // Make a fasta sequence
  // Use the first model
  const PDBModel & model = private_model[0];
  return model.to_fasta(get_id_code(), get_title());
}

vector<BioSequence> PDBEntry::to_coord_fasta() const {
  // Make a fasta sequence
  // Use the first model
  const PDBModel & model = private_model[0];
  return model.to_coord_fasta(get_id_code(), get_title());
}

vector<BioSequence> PDBEntry::to_seqres_fasta() const {
  // Make a fasta sequence
  // Use the first model
  const PDBModel & model = private_model[0];
  return model.to_seqres_fasta(get_id_code(), get_title());
}

void PDBEntry::clear() {
  private_model.clear();
  private_modres.clear();
  private_id_code = "????";
  private_classification = "";
  private_deposition_date = "";
  private_title = "";
}

istream & operator>>(istream & input_stream, PDBEntry & pdb)
{
  char line_buffer[1001];
  PDBModel model;
  pdb.clear();

  // Read one line at a time
  while (input_stream.good()) {
    input_stream.getline(line_buffer, 1000);
    string line = line_buffer;

    // ***********************************
    // *** Delegate Coordinate Entries ***
    // ***********************************
    // Quickly delegate frequent coordinate records to child object
    // (i.e. do this check first)
    if (line.substr(0,6) == "ATOM  ") {
      if (pdb.get_header_only()) break;
      model.add_line(line);
      continue;
    }
    if (line.substr(0,6) == "HETATM") {
      if (pdb.get_header_only()) break;
      model.add_line(line);
      continue;
    }
    
    // ***********************************
    // *** Parse Model/File Boundaries ***
    // ***********************************
    // Don't read past "END" cards
    if (line.substr(0,6) == "END   ") break;

    // Start a new model after ENDMDL card
    if (line.substr(0,6) == "ENDMDL") {
      if (pdb.get_header_only()) break;
      model.add_line(line);
      if (pdb.get_first_model_only()) break;
      // Only store the model if it has chains
      if (model.chain().size() > 0) pdb.add_model(model);
      model.clear();
      continue;
    }

    // ****************************
    // *** Parse Header records ***
    // ****************************
    if (line.substr(0,6) == "HEADER") {
      pdb.classification() = line.substr(10,40);
      pdb.deposition_date() = line.substr(50,9);
      pdb.id_code() = line.substr(62,4);
      continue;
    }
    // if (line.substr(0,6) == "TITLE ") {
	if (line.substr(0,6) == "COMPND") {
      pdb.title() += line.substr(10,60);
      // strip off trailing spaces
      string::iterator i = pdb.title().end();
      i --;
      while (*i == ' ') {
	i = pdb.title().erase(i);
	i --;
      }
      continue;
    }

    // Explicit residue sequence
    if (line.substr(0,6) == "SEQRES") {
      // int serial_number = atoi(line.substr(8, 2).c_str()); // not used?
      char chain_id = line.substr(11,1)[0];
      int num_residues = atoi(line.substr(13,4).c_str());

      map<char,Seqres> & s = pdb.seqres();
      if (s.find(chain_id) == s.end()) {
	Seqres new_seqres(chain_id);
	s[chain_id] = new_seqres;
      }
      Seqres & seqres = (*s.find(chain_id)).second;
      seqres.set_length(num_residues);

      // Loop over three-letter-code residue names
      int line_pos = 19;
      for (line_pos = 19; line_pos <= 68; line_pos += 4) {
	string residue_name = line.substr(line_pos, 3);
	if (residue_name.length() < 3) break;
	if (residue_name == "   ") continue;
	seqres.add_residue(residue_name);
      }
      continue;
    } 
    // COLUMNS        DATA TYPE       FIELD         DEFINITION
    // ---------------------------------------------------------------------------------
    //  1 -  6        Record name     "SEQRES"
    // 
    //  9 - 10        Integer         serNum        Serial number of the SEQRES record
    //                                              for the current chain.  Starts at 1
    //                                              and increments by one each line.
    //                                              Reset to 1 for each chain.
    // 
    // 12             Character       chainID       Chain identifier.  This may be any
    //                                              single legal character, including a
    //                                              blank which is used if there is
    //                                              only one chain.
    // 
    // 14 - 17        Integer         numRes        Number of residues in the chain.
    //                                              This value is repeated on every
    //                                              record.
    // 
    // 20 - 22        Residue name    resName       Residue name.
    // 
    // ...
    // 
    // 64 - 66        Residue name    resName       Residue name.
    // 
    // 68 - 70        Residue name    resName       Residue name.


    // FIXME parse DBREF cards
    // Sequence database cross references
    if (line.substr(0,6) == "DBREF ") {
      continue;
    }    

    // parse MODRES cards, and propagate info to PDBResidue objects
    // Parent normal residue types for modified residues
    if (line.substr(0,6) == "MODRES") {
      Modres modres(line);
      pdb.modres().push_back(modres);
      continue;
    }

    // *************************
    // *** Delegate the rest ***
    // *************************
    model.add_line(line); // Let PDBModel object handle other cases
  }
  // **********************************
  // *** Finished parsing all lines ***
  // **********************************

  // Always add the first model
  if ((pdb.get_model().size() < 1) ||
      (model.chain().size() > 0)) pdb.add_model(model);

  // ***********************************************
  // *** Attach MODRES info to selected residues ***
  // ***********************************************
  // Loop over modres records
  vector<Modres>::const_iterator modres_ptr;
  for (modres_ptr = pdb.modres().begin(); 
       modres_ptr != pdb.modres().end();
       modres_ptr ++) {
    // Loop over models
    vector<PDBModel>::iterator model_ptr;
    for (model_ptr = pdb.model().begin();
	 model_ptr != pdb.model().end();
	 model_ptr ++) {
      PDBChain & chain = model_ptr->chain(modres_ptr->get_chain_id());
      chain.set_modres(*modres_ptr);
    }
  }

  // Loop over seqres records
  map<char, Seqres>::const_iterator seqres_ptr;
  for (seqres_ptr = pdb.get_seqres().begin(); 
       seqres_ptr != pdb.get_seqres().end();
       seqres_ptr ++) {
    const Seqres & seqres = seqres_ptr->second;
    // if (seqres.error()) abort();
    // Loop over models
    vector<PDBModel>::iterator model_ptr;
    for (model_ptr = pdb.model().begin();
	 model_ptr != pdb.model().end();
	 model_ptr ++) {
      PDBChain & chain = model_ptr->chain(seqres.get_chain_id());
      chain.set_seqres(seqres);
    }
  }

  // Housekeeping of residues (check for redundancy)
  vector<PDBModel>::iterator model_ptr;
  for (model_ptr = pdb.model().begin();
       model_ptr != pdb.model().end();
       model_ptr ++) {
    model_ptr->check_residues();
  }

  return input_stream;
}

ostream & operator<<(ostream & os, const PDBEntry & pdb) {
  // Emit each model in turn
  vector<PDBModel>::const_iterator model_ptr;
  for (model_ptr = pdb.get_model().begin(); model_ptr != pdb.get_model().end(); ++ model_ptr) {
    os << *model_ptr;
  }
  // Finish up with the "END" card
  os << "END   " << endl;
  return os;
}
