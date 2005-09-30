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
// Revision 1.4  2002/09/13 23:28:08  bruns
// Added license header to all header files
//
// Revision 1.3  2002/05/15 19:14:02  bruns
// Added unary operator-
// Added length() prototype
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include "BaseVector.h"

class VectorND;
class Matrix3D;
class Vector3D;

class Vector3D : public BaseVector {
private:
  double private_element[3];
public:
  Vector3D() {private_element[0] = private_element[1] = private_element[2] = 0;}
  Vector3D(double x, double y, double z) {
    Vector3D & v = *this;
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }
  const double & get_x() const;
  const double & get_y() const;
  const double & get_z() const;

  double & set_x();
  double & set_y();
  double & set_z();
  double & set_x(double d) {set_x() = d; return set_x();}
  double & set_y(double d) {set_y() = d; return set_y();}
  double & set_z(double d) {set_z() = d; return set_z();}
  Vector3D & set(double x, double y, double z) {
    private_element[0] = x;
    private_element[1] = y;
    private_element[2] = z;
	return *this;
  }

  int get_size() const {return 3;}
  const double & operator[](const int index) const {return private_element[index];}
  double & operator[](const int index) {
    return private_element[index];
  }

  Vector3D operator+(const Vector3D & v2) const;
  Vector3D operator-(const Vector3D & v2) const;
  Vector3D operator*(double r) const;
  Vector3D operator-() const;

  Vector3D cross(const Vector3D & v2) const;


  // Numerical analysis
  Vector3D householder() const;  

  VectorND subvector(const int n1, const int n2) const;
  const double distance(const Vector3D & v2) const;
  const double distance_squared(const Vector3D & v2) const;
  Vector3D unit() const;
  double dot(const Vector3D & v2) const;
  double length() const;
};

#endif
