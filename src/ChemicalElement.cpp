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

