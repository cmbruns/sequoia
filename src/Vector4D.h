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
#ifndef _VECTOR4D_H_
#define _VECTOR4D_H_

#include "BaseVector.h"
class Matrix4D;
class VectorND;

class Vector4D : public BaseVector {
private:
  double private_element[4];
public:
  Vector4D() {private_element[0] = private_element[1] = private_element[2] = 0;}
  Vector4D & set(double w, double x, double y, double z) {
    private_element[0] = w;
    private_element[1] = x;
    private_element[2] = y;
    private_element[3] = z;
    return *this;
  }

  int get_size() const {return 4;}
  const double & operator[](const int index) const {return private_element[index];}
  double & operator[](const int index) {return private_element[index];}
  Vector4D operator+(const Vector4D & v2) const;
  Vector4D operator-(const Vector4D & v2) const;
  Vector4D operator*(double r) const;

  double dot(const Vector4D & v2) const;

  Vector4D unit() const;
  const double distance(const Vector4D & v2) const;
  VectorND subvector(const int n1, const int n2) const;
};

#endif
