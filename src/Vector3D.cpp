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
// Revision 1.4  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.3  2002/05/15 19:14:36  bruns
// Added routines for unary operator- and length()
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Vector3D.h"
#include <cmath>
#include "Matrix3D.h"
#include "VectorND.h"

#define SIGN(x) (((x) < 0) ? -1.0 : 1.0)
#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

const double & Vector3D::get_x() const {return private_element[0];}
const double & Vector3D::get_y() const {return private_element[1];}
const double & Vector3D::get_z() const {return private_element[2];}

double & Vector3D::set_x() {return private_element[0];}
double & Vector3D::set_y() {return private_element[1];}
double & Vector3D::set_z() {return private_element[2];}

const double Vector3D::distance(const Vector3D & v2) const {
    return sqrt(distance_squared(v2));
}

const double Vector3D::distance_squared(const Vector3D & v2) const {
    double sum = 0;
    double delta;
    delta = get_x() - v2.get_x();
    sum += delta * delta;
    delta = get_y() - v2.get_y();
    sum += delta * delta;
    delta = get_z() - v2.get_z();
    sum += delta * delta;
    return sum;
}

Vector3D Vector3D::operator+(const Vector3D & v2) const {
  Vector3D answer;
  int i;
  const Vector3D v1 = *this;

  answer = v1;
  for (i = 0; i < 3; ++i)
    answer[i] += v2[i];
  
  return answer;
}

Vector3D Vector3D::operator-(const Vector3D & v2) const {
  Vector3D answer;
  int i;
  const Vector3D & v1 = *this;

  answer = v1;
  for (i = 0; i < 3; ++i)
    answer[i] -= v2[i];
  
  return answer;
}

Vector3D Vector3D::operator-() const {
  Vector3D answer;
  int i;
  const Vector3D & v1 = *this;

  answer = v1;
  for (i = 0; i < 3; ++i)
    answer[i] = -answer[i];
  
  return answer;
}

Vector3D Vector3D::operator*(double r) const {
  Vector3D answer;
  int i;
  const Vector3D & v1 = *this;

  answer = v1;
  for (i = 0; i < 3; ++i)
    answer[i] = r * v1[i];
  
  return answer;
}

Vector3D Vector3D::unit() const {
  const Vector3D v = *this;
  double scale = v[0] * v[0] +
    v[1] * v[1] +
    v[2] * v[2];
  if (ABS(1.0 - scale) < EPSILON) return v;
  else {
    scale = sqrt(1.0 / scale);
    Vector3D answer = v * scale;
    return answer;
  }
}

double Vector3D::dot(const Vector3D & v2) const {
  double answer;
  int i;
  const Vector3D & v1 = *this;

  answer = 0;
  for (i = 0; i < 3; ++i)
    answer += v1[i] * v2[i];
  
  return answer;
}

Vector3D Vector3D::cross(const Vector3D & v2) const {
  Vector3D answer;
  const Vector3D & v1 = *this;

  answer[0] = v1[1] * v2[2] - v1[2] * v2[1];
  answer[1] = v1[2] * v2[0] - v1[0] * v2[2];
  answer[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return answer;
}

VectorND Vector3D::subvector(const int n1, const int n2) const {
	VectorND answer(n2 - n1 + 1);
	int i;
	// int j;
	const Vector3D & v = *this;

	int n = n2 - n1 + 1;

	for (i = 0; i < n; i++ )
		answer[i] = v[i + n1];

	return answer;
}

/* Householder vector - algorithm 5.1.1 in G&vL (p. 196) */
Vector3D Vector3D::householder() const {
  const Vector3D & x = *this;
  Vector3D v; /* answer */
  double mu; 
  double beta;
  int i;
  
  mu = sqrt(x.dot(x)); /* 2-norm of x */
  v = x;
  if (mu != 0) {
    beta = x[0] + (SIGN(x[0]) * mu);
    for (i = 1; i < 3; ++i)
      v[i] = v[i]/beta;
  }
  
  v[0] = 1.0;
  
  return v;
}
  
double Vector3D::length() const {
  return sqrt(this->dot(*this));
}

