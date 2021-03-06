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
// $Header$
// $Log$
// Revision 1.1  2004/06/04 19:34:47  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.7  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.6  2002/09/13 22:41:59  bruns
// Made get_seqres() have const qualifier
//
// Revision 1.5  2002/05/22 02:18:19  bruns
// Changed FastaSequence to BioSequence, to match new way
//
// Revision 1.4  2002/05/15 19:26:01  bruns
// Added get_first_residue method
//
// Revision 1.3  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
// Revision 1.2  2001/11/15 20:09:00  bruns
// Added cvs tags
//
#ifndef _PDBCHAIN_H_
#define _PDBCHAIN_H_

#include <vector>
#include <map>
#include <string>
#include "Seqres.h"
#include "PDBResidue.h"
#include "Hash3D.h"

class PDBAtom;
class BioSequence;
class Modres;

using namespace std;

// child of PDBModel, parent of PDBResidue
class PDBChain {
private:
  map<string, PDBResidue> private_residue;
  vector<string> private_residue_order;
  char private_chain_id;
  Seqres private_seqres;
  vector<Modres> private_modres;
  vector<Modres> private_used_modres; // MODRES cards that were actually assigned to coordinates
  string p_sequence_id;
  string p_title;
  double p_accessible_surface_area;
public:

	  double store_accessible_surface_area();
	  double compute_accessible_surface_area() const;
	  double get_accessible_surface_area() const {
		  if (p_accessible_surface_area < 0) throw AREA_NOT_YET_COMPUTED_EXCEPTION();
		  return p_accessible_surface_area;
	  }
	  
	  void assign_secondary_structure(); // TODO - assign by Kabsch Sander method

	  string & title() {return p_title;}
	  string & sequence_id() {return p_sequence_id;}
	  string get_title() const {return p_title;}
	  string get_sequence_id() const {return p_sequence_id;}
	  
	  bool is_protein(string residue_name);
	  const vector<string> & residue_order() const;
	  
	  void tweak_modres();
	  void set_modres(const Modres & modres);
	  PDBResidue * add_residue(const string & res_name,
							   const int & res_num,
							   const char & i_code);
	  void add_atom(const PDBAtom & atom);
	  void add_terminus(const string & res_name,
						const int & res_num,
						const char & i_code);
	  void set_seqres(const Seqres & seqres);
	  
	  PDBResidue * residue(const string & res_name,
						   const int & res_num,
						   const char & i_code);
	  PDBResidue & residue(unsigned int index) {
		  PDBResidue * res_ptr = residue_by_key(residue_order()[index]);
		  return *res_ptr;
	  }
	  PDBResidue * residue_by_key(const string & key);
	  PDBResidue & get_first_residue() {return *(residue_by_key(residue_order()[0]));}
	  
	  unsigned int residue_count() const {return residue_order().size();}
	  const PDBResidue & get_residue(unsigned int index) const {
		  return get_residue_by_key(residue_order()[index]);
	  }
	  const PDBResidue & get_residue(const string & res_name,
									 const int & res_num,
									 const char & i_code) const;
	  const PDBResidue & get_residue_by_key(const string & key) const;
	  const Seqres & get_seqres() const {return private_seqres;}
	  const char & get_chain_id() const;
	  
	  BioSequence to_fasta(const string & id, const string & title) const;
	  BioSequence to_coord_fasta(const string & id, const string & title) const;
	  BioSequence to_coord_fasta() const;
	  BioSequence to_seqres_fasta(const string & id, const string & title) const;
	  
	  void check_residues();
	  
	  PDBChain(char chain_id = ' ');
};

ostream & operator<<(ostream & os, const PDBChain & chain);

#endif
