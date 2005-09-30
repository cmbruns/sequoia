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
