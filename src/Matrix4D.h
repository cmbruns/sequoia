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
// Revision 1.1  2004/06/04 19:34:47  cmbruns
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
#ifndef _MATRIX4D_H_
#define _MATRIX4D_H_

// The assumption is that this class will usually be used for homogeneous
// rotation/translation matrices
// So we should try to keep them orthogonal

#include <iostream>
#include "BaseMatrix.h"
#include "Vector4D.h"

class Matrix3D;
class Vector3D;
class MatrixMN;

class Matrix4D : public BaseMatrix {
private:
  Vector4D private_row[4];
public:
  // Constructor
  Matrix4D() {
    // Default to identity matrix
    int i, j;
    for (i = 0; i < 4; i ++) 
      for (j = 0; j < 4; j ++) {
	if (i == j) private_row[i][j] = 1.0; // Diagonal
	else private_row[i][j] = 0.0; // Off-diagonal
      }
  }

  // Make homogeneous matrix from rotation matrix and translation
  Matrix4D & set(const Matrix3D & R, const Vector3D & t);

  const Vector4D & operator[](int index) const;
  Vector4D & operator[](int index); // Non-const version
  Matrix4D operator+(const Matrix4D & m2) const;
  Matrix4D operator-(const Matrix4D & m2) const;
  Matrix4D operator*(const Matrix4D & m2) const;
  Vector4D operator*(const Vector4D & v) const;
  Vector3D operator*(const Vector3D & v) const;

  double trace() const;
  Matrix4D scale(double scale) const;
  Matrix4D transpose() const;
  Matrix4D inverse_homog() const; // Inverse of homogeneous 3D transform matrix
  Matrix3D get_rotation() const;
  Vector3D get_translation() const;
  void set_rotation(const Matrix3D & R);
  void set_translation(const Vector3D & t);

  // FIXME - this could be moved to a common base class
  MatrixMN get_submatrix(const int m1, const int m2, const int n1, const int n2) const;
  Matrix4D & set_submatrix(const BaseMatrix & M2, int m, int n);
  int get_n_columns() const {return 4;}
  int get_n_rows() const {return 4;}

  // Numerical analysis
  Vector4D top_eigenvector() const;
  Vector4D sym_eigen(Matrix4D * Q = NULL) const;
  Matrix4D sym_tridiag(Matrix4D * Q = NULL) const;
  void row_rot(const int j, const int k, const Vector2D & cs);
  void col_rot(const int j, const int k, const Vector2D & cs);
};

Matrix4D eye4();
Matrix4D zeros4();

#endif
