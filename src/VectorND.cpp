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
// Revision 1.3  2002/09/13 22:17:36  bruns
// Removed redundant default function parameter value, to avoid warning from latest gcc
//
// Revision 1.2  2001/11/15 20:36:42  bruns
// Added cvs tags to [A-Z]*.cc and [A-Z]*.h
//
#include "VectorND.h"
#include "MatrixMN.h"

#define EPSILON FLT_EPSILON
#define ABS(x) (((x) < 0) ? (-(x)) : (x))
#define SIGN(x) (((x) < 0) ? -1.0 : 1.0)

VectorND::VectorND(int n) {
  int i;
  private_element.clear();
  for(i = 0; i < n; i ++)
    private_element.push_back(0);
}

int VectorND::get_size() const {return private_element.size();}

const double & VectorND::operator[](const int index) const {return private_element[index];}

double & VectorND::operator[](const int index) {return private_element[index];}

VectorND VectorND::operator+(const VectorND & v2) const {
  if (get_size() != v2.get_size()) abort();
  VectorND answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] += v2[i];
  return answer;
}

VectorND VectorND::operator-(const VectorND & v2) const {
  if (get_size() != v2.get_size()) abort();
  VectorND answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] -= v2[i];
  return answer;
}

VectorND VectorND::operator*(double r) const {
  VectorND answer = *this;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer[i] *= r;
  return answer;
}

double VectorND::dot(const VectorND & v2) const {
  const VectorND & v1 = *this;
  if (v1.get_size() != v2.get_size()) abort();
  double answer = 0.0;
  int i;
  for(i = 0; i < get_size(); i ++)
    answer += v1[i] * v2[i];
  return answer;
}

MatrixMN VectorND::outer_product(const VectorND & v2) const {
  const VectorND & v1 = *this;
  MatrixMN answer(v1.get_size(), v2.get_size());
  int i, j;
  
  for (i = 0; i < v1.get_size(); ++i) {
    for (j = 0; j < v2.get_size(); ++j)
      answer[i][j] = v1[i] * v2[j];
  }
  
  return answer;
}

VectorND VectorND::unit() const {
  const VectorND v = *this;
  double scale = v[0] * v[0] +
    v[1] * v[1] +
    v[2] * v[2];
  if (ABS(1.0 - scale) < EPSILON) return v;
  else {
    scale = sqrt(1.0 / scale);
    VectorND answer = v * scale;
    return answer;
  }
}

const double VectorND::distance(const VectorND & v2) const {
  const VectorND & v1 = *this;
  VectorND delta = v2 - v1;
  return sqrt(delta.dot(delta));
}

/* Householder vector - algorithm 5.1.1 in G&vL (p. 196) */
VectorND VectorND::householder() const {
  const VectorND & x = *this;
  VectorND v = *this; /* answer */
  double mu; 
  double beta;
  int i;

  mu = sqrt(x.dot(x)); /* 2-norm of x */
  if (mu != 0) {
    beta = x[0] + (SIGN(x[0]) * mu);
    for (i = 1; i < x.get_size(); ++i)
      v[i] = v[i]/beta;
  }

  v[0] = 1.0;

  return v;
}

