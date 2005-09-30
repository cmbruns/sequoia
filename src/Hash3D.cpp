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
// Revision 1.1  2004/06/04 19:34:48  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.3  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.2  2002/05/23 18:04:44  bruns
// Added CVS tags
//
#include "Hash3D.h"
#include "PDBAtom.h"
#include "PDBResidue.h"
#include "PDBChain.h"
#include "ChemicalElement.h"

Vector3D PDBAtomHash3DContent::get_position() const {
	return p_atom->get_coordinate();
}

void Hash3D::add_pdb_residue(const PDBResidue & residue) {
	for (unsigned int a = 0; a < residue.get_atom().size(); a++) {
		// TODO - restrict to heavy protein atoms
		const PDBAtom & atom = residue.get_atom()[a];
		if (atom.get_record_name() != "ATOM  ") continue; // Protein only
		if (atom.get_element() == & HYDROGEN) continue; // No hydrogen
		if (atom.get_element() == NULL) continue; // No unknown
		add_pdb_atom(atom);
	}
}

void Hash3D::add_pdb_chain(const PDBChain & chain) {
	for (unsigned int r = 0; r < chain.residue_count(); r++) {
		// TODO - restrict to protein residues
		const PDBResidue & residue = chain.get_residue(r);
		if (! residue.is_protein()) continue; // protein only
		// Is the first residue always "redundant"?
		// if (residue.is_redundant()) continue; // just one version of each residue
		add_pdb_residue(residue);
	}
}


