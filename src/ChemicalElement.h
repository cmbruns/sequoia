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
