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
