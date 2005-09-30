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
