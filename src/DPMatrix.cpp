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

#include "DPMatrix.h"

ostream & DPMatrix::print_debug(ostream & os, unsigned int indent_size) const {
	const DPMatrix & mat = *this;
	string indent = "";
	for (unsigned int i=0;i<indent_size;i++) {indent += " ";}
	
	for (unsigned int i = 0; i < mat.size(); i++) {
		for (unsigned int j = 0; j < mat[i].size(); j++) {
			os << indent << "Cell(" << i << ", " << j << "):" << endl;
			mat[i][j].print_debug(os, indent_size + 2);
			os << endl;
		}
	}
	
	return os;
}

