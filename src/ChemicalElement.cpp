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
// Revision 1.1  2004/06/04 19:34:48  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//

#include "ChemicalElement.h"

#include <algorithm>
#include <cctype>
#include "Exceptions.h"

const ChemicalElement CARBON("C", 1.90);
const ChemicalElement HYDROGEN("H", 0.00);
const ChemicalElement NITROGEN("N", 1.54);
const ChemicalElement OXYGEN("O", 1.48);
const ChemicalElement SULFUR("S", 1.65);
const ChemicalElement UNKNOWN_ELEMENT("Xx", 0.00);

bool ChemicalElement::operator==(const ChemicalElement & element2) const {
	if ((symbol.length() > 0) && 
		(element2.symbol.length() > 0) &&
		(symbol == element2.symbol)) 
		return true;
	return false;
}

char myupper(const char & c) {return toupper(c);}

const ChemicalElement * find_element(string symbol) {
	// TODO - massage element symbol
	string s = symbol;
	if (s.length() == 4) // Probably a PDB atom name
		s = s.substr(0,2); // Element should be in first 2 characters
	if (s.length() > 2) throw STRANGE_ELEMENT_NAME_EXCEPTION();
	
	// Change everything into upper case
	transform(s.begin(), s.end(), s.begin(), myupper);
	
	// Keep only the letters
	int element_start = s.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	if (element_start == string::npos) throw STRANGE_ELEMENT_NAME_EXCEPTION();
	int element_end = s.find_last_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	s = s.substr(element_start, element_end - element_start + 1);
	
	if (s == "C") return & CARBON;
	if (s == "O") return & OXYGEN;
	if (s == "H") return & HYDROGEN;
	if (s == "N") return & NITROGEN;
	if (s == "S") return & SULFUR;
	return NULL; // Unknown Element
}

