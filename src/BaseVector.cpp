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
// Revision 1.1  2004/06/04 19:34:46  cmbruns
// Imported structure related sources from archive on baxter
// Debugged simple conversion of structures to sequences.
// Implemented computation of solvent accessible surface areas
// Created target residue_area, for output of residue solvent accessible surfaces areas
// Updated GPL headers
//
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2001/12/14 19:21:55  bruns
// Changes to compile with gcc 3.0.2
//  - "using namespace std;"  to use standard classes
//  - "ios::fmtflags" instead of "fmtflags"
//
// Revision 1.2  2001/11/15 20:19:55  bruns
// Added cvs tags up through FastaSequence.h
//
#include "BaseVector.h"
#include "MatrixMN.h"

MatrixMN BaseVector::to_row() const {
  MatrixMN answer(1, get_size());
  const BaseVector & v = *this;
  int i;
  for (i = 0; i < get_size(); i++)
    answer[0][i] = v[i];
    return answer;
}
MatrixMN BaseVector::to_column() const {
  MatrixMN answer(get_size(), 1);
  const BaseVector & v = *this;
  int i;
  for (i = 0; i < get_size(); i++)
    answer[i][0] = v[i];
  return answer;
}
	
ostream & operator<<(ostream & os, const BaseVector & v) {
  // Remember initial format flags settings
  ios::fmtflags old_format = os.flags(); // new iostream way
  // long old_format = os.flags(); // old iostream way

  int i;
  os << "(";
  for (i = 0; i < v.get_size(); ++i) {
    os.setf(ios::fixed); // Makes precision be what I want?
    os.precision(5); os.width(9);
    os << v[i];
    if (i < (v.get_size() - 1)) os << "   ";
  }
  os << ")";

  // Restore initial settings
  os.flags(old_format);

  return os;
}
