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
