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
// $Header$
// $Log$
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//

#ifndef __CHEMICAL_ELEMENT_H__
#define __CHEMICAL_ELEMENT_H__

#include <string>

using namespace std;

class ChemicalElement {
public:
	string symbol;
	// double mass; // Not used yet
	double vdw_radius; // in Angstroms, for use in solvent accessibility
	bool operator==(const ChemicalElement & element2) const;
	ChemicalElement(const string & symbol_string, double radius) : 
		symbol(symbol_string), vdw_radius(radius)
	{
	}
};

const ChemicalElement * find_element(string symbol);

extern const ChemicalElement CARBON;
extern const ChemicalElement HYDROGEN;
extern const ChemicalElement NITROGEN;
extern const ChemicalElement OXYGEN;
extern const ChemicalElement SULFUR;
extern const ChemicalElement UNKNOWN_ELEMENT;

#endif
