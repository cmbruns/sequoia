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
// Revision 1.3  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _MATRIXMN_H_
#define _MATRIXMN_H_

class MatrixMN;

#include <vector>
#include "BaseMatrix.h"
#include "VectorND.h"

class Vector2D;

class MatrixMN : public BaseMatrix {
private:
  vector<VectorND> private_element;
protected:
public:
  int get_n_rows() const;
  int get_n_columns() const;

  VectorND & operator[](int index);
  const VectorND & operator[](int index) const;
  VectorND operator*(const VectorND & v) const;
  MatrixMN operator+(const MatrixMN & m2) const;
  MatrixMN transpose() const;

  // Numerical analysis
  MatrixMN sym_QRstep(int p, BaseMatrix * Q = NULL) const;

  void row_rot(const int i, const int k, const Vector2D & cs);
  void col_rot(const int i, const int k, const Vector2D & cs);

  MatrixMN(){}
  MatrixMN(int rows, int cols) {
	  int i;
	  VectorND v(cols);
	  for (i = 0; i < rows; i ++) {
		  private_element.push_back(v);
	  }
  }
  virtual ~MatrixMN() {}
};

#endif
