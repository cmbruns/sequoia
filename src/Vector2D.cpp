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
// Revision 1.3  2002/09/14 00:02:51  bruns
// Added license header to most .cc files
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "Vector2D.h"

#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))

Vector2D::Vector2D() {
  int i;
  for(i = 0; i < 2; i ++)
    private_element[i] = 0;
}
int Vector2D::get_size() const {return 2;}

const double & Vector2D::operator[](const int index) const {return private_element[index];}
double & Vector2D::operator[](const int index) {return private_element[index];}
Vector2D Vector2D::operator+(const Vector2D & v2) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] += v2[i];
  return answer;
}
Vector2D Vector2D::operator-(const Vector2D & v2) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] -= v2[i];
  return answer;
}
Vector2D Vector2D::operator*(double r) const {
  Vector2D answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] *= r;
  return answer;
}
double Vector2D::dot(const Vector2D & v2) const {
  const Vector2D & v1 = *this;
  double answer = 0.0;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer += v1[i] * v2[i];
  return answer;
}

Vector2D Vector2D::unit() const {
  const Vector2D v = *this;
  double scale = v[0] * v[0] +
    v[1] * v[1] +
    v[2] * v[2];
  if (ABS(1.0 - scale) < EPSILON) return v;
  else {
    scale = sqrt(1.0 / scale);
    Vector2D answer = v * scale;
    return answer;
  }
}

const double Vector2D::distance(const Vector2D & v2) const {
  const Vector2D & v1 = *this;
  Vector2D delta = v2 - v1;
  return sqrt(delta.dot(delta));
}

Vector2D givens(const double a, const double b) {
  Vector2D answer;
  double tau, c, s;

  if (b == 0) {
    c = 1.0;
    s = 0.0;
  }
  else {
    if (ABS(b) > ABS(a)) {
      tau = -a / b;
      s = 1.0 / sqrt(1 + tau*tau);
      c = s * tau;
    }
    else {
      tau = -b / a;
      c = 1.0 / sqrt(1.0 + tau*tau);
      s = c * tau;
    }
  }
  
  answer[0] = c;
  answer[1] = s;
  return answer;
}

