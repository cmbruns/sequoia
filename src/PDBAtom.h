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
// Revision 1.4  2002/09/13 22:43:01  bruns
// Added copyright header
//
// Revision 1.3  2002/05/15 19:27:04  bruns
// Added get_coordinate() accessor, and get_dihedral_angle subroutine
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//

#ifndef _PDBATOM_H_
#define _PDBATOM_H_

#include <iostream>
#include <string>
#include <cmath>
#include "Vector3D.h"
#include "Exceptions.h"

using namespace std;

class ChemicalElement;
class Hash3D;

// child of PDBResidue, parent of no one
class PDBAtom {
private:
  // Standard fields
  string   private_record_name;
  int      private_serial_number;
  string   private_atom_name;
  char     private_alt_loc;
  string   private_residue_name;
  char     private_chain_id;
  int      private_residue_number;
  char     private_insertion_code;
  Vector3D private_coordinate;
  float    private_occupancy;
  float    private_temperature_factor;
  string   private_segment_identifier;
  string   private_element_symbol;
  string   private_charge;
  
  // Derived fields
  const ChemicalElement * p_element;
  double p_accessible_surface_area;
  
  void set_element();
  
public:
  string   get_record_name() const;
  int      get_serial_number() const;
  string   get_atom_name() const;
  char     get_alt_loc() const;
  string   get_residue_name() const;
  char     get_chain_id() const;
  int      get_residue_number() const;
  char     get_insertion_code() const;
  Vector3D get_coordinate() const;
  float    get_occupancy() const;
  float    get_temperature_factor() const;
  string   get_segment_identifier() const;
  string   get_element_symbol() const;
  string   get_charge() const;

  const ChemicalElement * get_element() const {return p_element;}
  double get_accessible_surface_area() const {
	  if (p_accessible_surface_area < 0) throw AREA_NOT_YET_COMPUTED_EXCEPTION();
	  return p_accessible_surface_area;
  }
  
  // Direct accessor
  Vector3D & coordinate() {
    return private_coordinate;
  }

  double distance(const PDBAtom & atom2) const;
  double compute_accessible_surface_area(const Hash3D & hash3d) const;
  double store_accessible_surface_area(const Hash3D & hash3d);

  PDBAtom() : p_element(NULL), p_accessible_surface_area(-1) {}
  PDBAtom(const string & s);
};

ostream & operator<<(ostream & os, const PDBAtom & atom);

double get_dihedral_angle(const PDBAtom & atom1,
		   const PDBAtom & atom2,
		   const PDBAtom & atom3,
		   const PDBAtom & atom4);
#endif
