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
// Revision 1.4  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.3  2002/05/15 19:28:11  bruns
// Added axis_angle prototype
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _MATRIX3D_H_
#define _MATRIX3D_H_

// The assumption is that this class will usually be used for rotation matrices
// So we should try to keep them orthogonal

#include <iostream>
#include "BaseMatrix.h"
#include "Vector3D.h"

class Quaternion;
class MatrixMN;

class Matrix3D : public BaseMatrix {
private:
  Vector3D private_row[3];
public:
  // Constructor
  Matrix3D();
  Matrix3D(double e00, double e01, double e02,
	   double e10, double e11, double e12,
	   double e20, double e21, double e22);
  Matrix3D & set(const Vector3D & v1, 
		 const Vector3D & v2,
		 const Vector3D & v3);
  int get_n_columns() const {return 3;}
  int get_n_rows() const {return 3;}

  const Vector3D & operator[](int index) const;
  Vector3D & operator[](int index); // Non-const version

  Matrix3D operator+(const Matrix3D & m2) const;
  Matrix3D operator-(const Matrix3D & m2) const;
  Matrix3D operator*(const Matrix3D & m2) const;
  Vector3D operator*(const Vector3D & v) const;
  Matrix3D operator*(const double scale) const;

  double trace() const;
  Matrix3D scale(double scale) const;
  Matrix3D transpose() const;
  Quaternion quaternion() const;
  
  MatrixMN get_submatrix(const int m1, const int m2, const int n1, const int n2) const;
  Matrix3D & Matrix3D::set_submatrix(const MatrixMN & M2, int m, int n);

  // Numerical analysis
  Vector3D top_eigenvector() const;
  Vector3D sym_eigen(Matrix3D * Q = NULL) const;
  void row_rot(const int i, const int k, const Vector2D & cs);
  void col_rot(const int i, const int k, const Vector2D & cs);
  Matrix3D sym_tridiag(Matrix3D * Q = NULL) const;
};

Matrix3D zeros3();
Matrix3D eye3();
Matrix3D axis_angle(const Vector3D & axis, const double angle);

#endif
