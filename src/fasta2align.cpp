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
// Revision 1.5  2004/06/04 19:32:24  cmbruns
// Updated GPL header
//
// Changed for initial testing of structure accessible surface area computation.  This has been extended in residue_area.cpp
//


#include <iostream>
#include <fstream>
#include "AlignmentMethod.h"
#include "BioSequence.h"
#include "Exceptions.h"

using namespace std;

int main (int argc, char * const argv[]) {

	// TODO - make this a command line parameter
	const char * = "pdb1coa.ent";
	
	PDBEntry test_pdb;
	test_pdb.load_pdb_file(pdb_file_name);
	
	PDBModel & pdb_model = test_pdb.get_first_model();
	for (char chain_id = pdb_model.chain_order().begin();
		 chain_id != pdb_modle.chain_order().end();
		 chain_id ++) {
		
		PDBChain & chain = pdb_model.chain(chain_id);

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
	}
		
	return 0;
	
	SequenceAlignment a1;
	a1.load_pdb_file("1FND.pdb");
	a1.set_weight(1.0);
	
	cout << a1 << endl;
	
	// Read second sequence file
	SequenceAlignment a2;
	a2.load_fasta_file("pdr.fasta");
	a2.set_weight(0.7);
	
	cout << a1.align(a2);
	
    return 0;
}
