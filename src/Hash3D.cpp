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


