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
// Revision 1.8  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.7  2002/05/22 02:13:40  bruns
// Added set_residue_number() method, now required by BaseResidue class
//
// Revision 1.6  2002/05/15 19:20:05  bruns
// Changed get_atom() to return an atom, not a pointer (uses modern exception
//  handling now)
// Added incomplete get_torsion routine, which presently only works with leucine
// side chains
//
// Revision 1.5  2001/11/19 20:47:01  bruns
// Added required new_clone() function ("virtual constructor")
//
// Revision 1.4  2001/11/16 00:16:19  bruns
// Added is_gap function, required by base class
//
// Revision 1.3  2001/11/15 21:36:30  bruns
// Made PDBResidue a descendant of new BaseResidue class
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _PDBRESIDUE_H_
#define _PDBRESIDUE_H_

#include <vector>
#include "Exceptions.h"
#include "Residue.h"
#include "PDBAtom.h"
#include "Modres.h"

class Hash3D;

string generate_unique_residue_key(string residue_name, 
				   int residue_number,
				   char insertion_code);
bool three_letter_is_nucleotide(string residue_name);
bool three_letter_is_protein(string residue_name);
char three_letter_to_one_letter_code(string residue_name);

// child of PDBChain, parent of PDBAtom
class PDBResidue : public Residue {
private:
  vector<PDBAtom> private_atom;
  string private_residue_name;
  string private_standard_residue_name; // As indicated on MODRES cards
  int private_residue_number;
  char private_insertion_code;
  string private_unique_key;
  char private_chain_id;

  bool private_is_hetatm_only;
  bool private_has_modres;
  bool private_is_terminus;
  bool private_is_protein;
  bool private_has_altloc; // non-blank alt_loc character for any atoms?
  bool private_is_redundant; // Is the the "same" as the previous residue, for sequence purposes?
  char private_one_letter_code;
  double p_accessible_surface_area;

public:
  const char & one_letter_code() const;
  bool is_gap() const {return false;}
  PDBResidue * new_clone() const {return new PDBResidue(get_chain_id(),
                                                        get_residue_name(),
                                                        get_residue_number(),
                                                        get_insertion_code());}

  void set_one_letter_code();
  void set_residue_number(int n) {private_residue_number = n;}

  // "get_" methods are const, plain ones are not
  PDBAtom & atom(string atom_name);
  vector<PDBAtom> & atom();
  const vector<PDBAtom> & get_atom() const;
  const PDBAtom get_atom(string atom_name) const;
  const string & get_residue_name() const {return private_residue_name;}
  const string & get_standard_residue_name() const {return private_standard_residue_name;}
  int get_residue_number() const {return private_residue_number;}
  char get_insertion_code() const {return private_insertion_code;}
  char get_chain_id() const {return private_chain_id;}

  int set_modres(const Modres & modres);

  bool is_protein() const;
  bool is_hetatm_only() const;
  bool has_modres() const;
  bool is_terminus() const;
  bool has_altloc() const;
  bool is_redundant() const;

  void set_terminus();
  void add_atom(const PDBAtom & atom);
  void check_redundancy(const PDBResidue * previous_residue);
  double compute_accessible_surface_area(const Hash3D & hash3d) const;
  double store_accessible_surface_area(const Hash3D & hash3d);
  double get_accessible_surface_area() const {
	  if (p_accessible_surface_area < 0) throw AREA_NOT_YET_COMPUTED_EXCEPTION();
	  return p_accessible_surface_area;
  }
  
  double get_torsion(const string & name) const {
    try {
      if ((name == "CHI1") || 
	  (name == "Chi1") ||
	  (name == "chi1") ||
	  (name == "x1") ||
	  (name == "X1")
	  ) {
	return get_dihedral_angle(get_atom(" N  "),
				  get_atom(" CA "),
				  get_atom(" CB "),
				  get_atom(" CG ")
				  );
      }
      if ((name == "CHI2") || 
	  (name == "Chi2") ||
	  (name == "chi2") ||
	  (name == "x2") ||
	  (name == "X2")
	  ) {
	return get_dihedral_angle(
				  get_atom(" CA "),
				  get_atom(" CB "),
				  get_atom(" CG "),
				  get_atom(" CD1")
				  );
      }
    } catch (no_such_atom_error e) {
      // cerr << "No such atom " << e.get_atom_name() << endl;
    }
    throw no_such_torsion_error();
  }

  void set_torsion(const string & torsion_name, const double angle);

  PDBResidue(char chain_id = ' ', string residue_name = "UNK", int residue_number = -1, char insertion_code = ' ');
  virtual ~PDBResidue() {}
};

ostream & operator<<(ostream & os, const PDBResidue & residue);

#endif
