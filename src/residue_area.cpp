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


#include <iostream>
#include "Exceptions.h"
#include "PDBEntry.h"

using namespace std;

int main (int argc, char * const argv[]) {
	
	// read file from command line parameter
	PDBEntry pdb_entry;
	if (argc > 1) {
		const char * pdb_file_name = argv[1];	
		cout << "Reading structure from file " << pdb_file_name << endl;
		pdb_entry.load_pdb_file(pdb_file_name);
	}
	else {
		cout << "Reading structure from standard input." << endl;
		cin >> pdb_entry;
	}
	
	PDBModel & pdb_model = pdb_entry.get_first_model();
	for (vector<char>::const_iterator chain_i = pdb_model.chain_order().begin();
		 chain_i != pdb_model.chain_order().end();
		 chain_i ++) {
		const char chain_id = *chain_i;
		
		PDBChain & chain = pdb_model.chain(chain_id);
		cout << "Chain ID = \"" << chain_id << "\"" << endl;
		
		chain.store_accessible_surface_area();
		for (unsigned int r = 0; r < chain.residue_count(); r++) {
			const PDBResidue & residue = chain.get_residue(r);
			
			double residue_area = 0.0;
			try {residue_area = residue.get_accessible_surface_area();}
			catch (AREA_NOT_YET_COMPUTED_EXCEPTION e) {continue;}
			
			cout << residue.one_letter_code() << " ";
			cout << residue.get_residue_number() << " ";
			cout << residue_area << endl;
		}
		cout << "Total chain protein accessible area = " << chain.get_accessible_surface_area() << endl;
		
	}
	
	return 0;
}
