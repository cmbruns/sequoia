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


#ifndef __DPMATRIX_H__
#define __DPMATRIX_H__

#include <vector>
#include "ConservidueAlignment.h"

class DPRow : public vector<ConservidueAlignment> {
public:
	DPRow(unsigned int size, const ConservidueAlignment & cell) : vector<ConservidueAlignment>(size, cell) {}
};

class DPMatrix : public vector<DPRow> {
public:
	DPMatrix(unsigned int size, const DPRow & row) : vector<DPRow>(size, row) {}
	ostream & print_debug(ostream & os, unsigned int indent_size) const;
};

#endif