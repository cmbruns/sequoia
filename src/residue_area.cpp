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
